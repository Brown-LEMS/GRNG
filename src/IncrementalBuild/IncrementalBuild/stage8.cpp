/*
Copyright 2022, Brown University, Providence, RI.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose other than its incorporation into a
commercial product or service is hereby granted without fee, provided
that the above copyright notice appear in all copies and that both
that copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
// Cole Foster
// 2022-06-15
#include "incremental-build.hpp"
#include "omp.h"

namespace IncrementalBuild_Stage8 {
void collectNeighborsChildren(QueryStructIncrementalBuild& queryStruct, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                              tsl::sparse_set<unsigned int>& candidateNeighborList);

bool checkLinkForInterference(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex,
                              std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

void collectCandidateCoarseIndices(unsigned int const queryIndex, int const layerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                   tsl::sparse_set<unsigned int> const& activeInterferingPivotsList, RankOrderedList& rankedCandidateList);

bool dmaxCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance, float const activeChildLinkDistance,
               int const layerIndex, int const finerLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
               RankOrderedList const& candidateInterferingPivots, tsl::sparse_set<unsigned int>& activeInterferingList);

void collectCandidateFineIndices(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                 float const activeChildLinkDistance, int const finerLayerIndex,
                                 std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                 tsl::sparse_set<unsigned int> const& activeCoarseInterferingPivotsList,
                                 tsl::sparse_set<unsigned int>& activeFineInterferingPivotsList);

bool assistedChildCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                        float const activeChildLinkDistance, unsigned int const coarseInterferingPivot, unsigned int const fineInterferingPivot,
                        std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

bool bruteForceInterferenceCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                 float const activeChildLinkDistance, tsl::sparse_set<unsigned int> const& interferingPivotsList,
                                 std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

};  // namespace IncrementalBuild_Stage8

void IncrementalBuild_Stage8::collectNeighborsChildren(QueryStructIncrementalBuild& queryStruct,
                                                       std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                       tsl::sparse_set<unsigned int>& candidateNeighborList) {
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    int const finerLayerIndex = pivotLayerIndex + 1;

    // collect children of all neighbors as potential virtual neighbors
    tsl::sparse_set<unsigned int> activeChildrenList{};
    tsl::sparse_set<unsigned int> const& pivotLayerNeighbors = queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex];
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotLayerNeighbors.begin(); it1 != pivotLayerNeighbors.end(); it1++) {
        unsigned int const neighborIndex = (*it1);

        tsl::sparse_set<unsigned int> const& pivotNeighborChildren =
            (*pivotLayers)[pivotLayerIndex].get_descendantPivotIndices(neighborIndex, finerLayerIndex);
        activeChildrenList.insert(pivotNeighborChildren.begin(), pivotNeighborChildren.end());
    }

    // make sure all parents of neighbors of query
    candidateNeighborList.clear();
    for (it1 = activeChildrenList.begin(); it1 != activeChildrenList.end(); it1++) {
        unsigned int const childIndex = (*it1);
        bool flag_addChild = true;

        // ensure query is neighbors to all parents of child
        tsl::sparse_set<unsigned int> const& childParentList = (*pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(childIndex, pivotLayerIndex);
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = childParentList.begin(); it2 != childParentList.end(); it2++) {
            unsigned int const childParentIndex = (*it2);

            if (pivotLayerNeighbors.find(childParentIndex) == pivotLayerNeighbors.end()) {
                flag_addChild = false;
                break;
            }
        }

        if (flag_addChild) {
            candidateNeighborList.insert(childIndex);
        }
    }
    return;
}

void IncrementalBuild::stage8(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int>& finerLayerPivotNeighborsList = queryStruct.finerCoarseNeighborsList;

    // collect activeChildren as children of neighbors
    tsl::sparse_set<unsigned int> activeChildrenList{};
    IncrementalBuild_Stage8::collectNeighborsChildren(queryStruct, _pivotLayers, activeChildrenList);

    // top-down, exhaustive check for interference of link(query,activeChild)
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeChildrenList.begin(); it1 != activeChildrenList.end(); it1++) {
        unsigned int const activeChildIndex = (*it1);

        bool flag_interference = IncrementalBuild_Stage8::checkLinkForInterference(queryStruct, activeChildIndex, _pivotLayers);
        if (!flag_interference) {
            finerLayerPivotNeighborsList.insert(activeChildIndex);
        }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[8] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[8] += (dEnd_stage - dStart_stage);
    return;
}

/**
 * @brief Exhaustively checks potential link(query,activeChildIndex) for interference in top-down manner
 *
 * @param queryStruct
 * @param activeChildIndex
 * @param finerLayerIndex
 * @param pivotLayers
 * @return true
 * @return false
 */
bool IncrementalBuild_Stage8::checkLinkForInterference(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex,
                                                       std::shared_ptr<PivotLayerHierarchy> const& pivotLayers) {
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    float const activeChildRadius = (*pivotLayers)[pivotLayerIndex + 1].get_pivotRadius(activeChildIndex);
    float const distance12 = pivotLayers->getDistance(queryIndex, activeChildIndex);
    float queryLinkDistance = pivotLayers->_luneRadius(distance12, queryRadius, activeChildRadius);
    float activeChildLinkDistance = pivotLayers->_luneRadius(distance12, activeChildRadius, queryRadius);

    // no lune, return no intererence
    if (queryLinkDistance <= 0 || activeChildLinkDistance <= 0) {
        return false;
    }

    bool flag_interference = false;

    // //------------------------------------------------------------------------------
    // // Work top-down to collect all coarse pivots whose children could interfere
    // //------------------------------------------------------------------------------
    tsl::sparse_set<unsigned int> interferingPivotsList{};
    if (pivotLayerIndex == 0) {
        tsl::sparse_set<unsigned int> const& indices = *(*pivotLayers)[pivotLayerIndex].get_pivotIndices_ptr();
        interferingPivotsList = indices;
    } else {
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>> interferingPivotsMap;
        for (int layerIndex = 0; layerIndex < pivotLayerIndex; layerIndex++) {
            interferingPivotsMap[layerIndex] = tsl::sparse_set<unsigned int>{};

            // created a ranked-list of candidate interfering pivots
            RankOrderedList candidateInterferingPivots;
            if (layerIndex == 0) {
                candidateInterferingPivots = queryStruct.rankOrderedTopLayerPivotsList;
            } else {
                collectCandidateCoarseIndices(queryIndex, layerIndex, pivotLayers, interferingPivotsMap[layerIndex - 1], candidateInterferingPivots);
            }

            // create list of actively interfering pivots from the candidate interfering using dmax
            // return true if any pivots found to be interfering
            flag_interference = dmaxCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, layerIndex, pivotLayerIndex,
                                          pivotLayers, candidateInterferingPivots, interferingPivotsMap[layerIndex]);

            if (flag_interference == true) {
                return true;  // interference, this link is no longer to be considered!
            }
        }

        //------------------------------------------------------------------------------
        // Collect potentially interfering fine pivots from children of interfering coarse pivots
        //------------------------------------------------------------------------------
        collectCandidateFineIndices(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, pivotLayerIndex, pivotLayers,
                                    interferingPivotsMap[pivotLayerIndex - 1], interferingPivotsList);
    }

    //------------------------------------------------------------------------------
    // Finally, brute force check for interference
    //------------------------------------------------------------------------------
    flag_interference =
        bruteForceInterferenceCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, interferingPivotsList, pivotLayers);

    return flag_interference;
}

/**
 * @brief Collect potentially interfering fine pivots as children of coarse potentially interfering
 * collectCandidateCoarseIndices(queryIndex, layerIndex, pivotLayers, interferingPivotsMap[layerIndex - 1], candidateInterferingPivots);
 * @param queryIndex
 * @param layerIndex
 * @param pivotLayers
 * @param activeInterferingPivotsList
 * @param rankedCandidateList
 */
void IncrementalBuild_Stage8::collectCandidateCoarseIndices(unsigned int const queryIndex, int const layerIndex,
                                                            std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                            tsl::sparse_set<unsigned int> const& activeInterferingPivotsList,
                                                            RankOrderedList& rankedCandidateList) {
    rankedCandidateList.clear();

    // iterate through each active pivot in layerIndex-1
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeInterferingPivotsList.begin(); it1 != activeInterferingPivotsList.end(); it1++) {
        unsigned int const activeIndex = (*it1);
        tsl::sparse_set<unsigned int> const& activeIndexChildrenList =
            (*pivotLayers)[layerIndex - 1].get_descendantPivotIndices(activeIndex, layerIndex);

        // consider all children of activeIndex as potential candidate
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activeIndexChildrenList.begin(); it2 != activeIndexChildrenList.end(); it2++) {
            unsigned int const childIndex = (*it2);
            bool flag_add = true;

            // check that all parents of childIndex are in activeInterferingPivotsList
            tsl::sparse_set<unsigned int> const& childIndexParentsList =
                (*pivotLayers)[layerIndex].get_ancestorPivotIndices(childIndex, layerIndex - 1);
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = childIndexParentsList.begin(); it3 != childIndexParentsList.end(); it3++) {
                unsigned int const childParentIndex = (*it3);

                if (activeInterferingPivotsList.find(childParentIndex) == activeInterferingPivotsList.end()) {
                    flag_add = false;
                    break;
                }
            }

            // child becomes potentially interfering pivot
            if (flag_add) {
                float const distance = (*pivotLayers).getDistance(queryIndex, childIndex);
                rankedCandidateList.add(childIndex, distance);
            }
        }
    }
    rankedCandidateList.sort();
}

/**
 * @brief Given list of potentially interfering pivots, use dmax to tell if any children within their domain may interfere with link.
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param layerIndex
 * @param finerLayerIndex
 * @param pivotLayers
 * @param candidateInterferingPivots
 * @param activeInterferingList
 * @return true
 * @return false
 */
bool IncrementalBuild_Stage8::dmaxCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                        float const activeChildLinkDistance, int const layerIndex, int const pivotLayerIndex,
                                        std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, RankOrderedList const& candidateInterferingPivots,
                                        tsl::sparse_set<unsigned int>& activeInterferingList) {
    activeInterferingList.clear();

    // create list of actively interfering pivots
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = candidateInterferingPivots.list.begin(); it1 != candidateInterferingPivots.list.end(); it1++) {
        unsigned int const candidateInterferingPivotIndex = (*it1).index;
        float const distance13 = (*it1).distance;
        bool flag_add = true;

        // // check lune side1 against maxDmax, the max of all possible parent-child distances
        float const dmaxMaxLayer = (*pivotLayers)[layerIndex].get_dmaxMaximumOfMaxChildDistance(pivotLayerIndex);
        if ((distance13 - dmaxMaxLayer) >= queryLinkDistance) {
            break;  // break since ranked list
        }

        // check lune side1 against dmax
        float const dmaxCandidateInterferingIndex =
            (*pivotLayers)[layerIndex].get_dmaxMaxChildDistance(candidateInterferingPivotIndex, pivotLayerIndex);
        if ((distance13 - dmaxCandidateInterferingIndex) >= queryLinkDistance) {
            flag_add = false;
            continue;  // no children of this one can interfere, next one!
        }

        // if distance available, lune side 2 with dmax and brute force check
        std::pair<bool, float> distance23_pair = (*pivotLayers).getDistanceIfAvailable(activeChildIndex, candidateInterferingPivotIndex);
        if (distance23_pair.first) {
            
            // brute force check for pivot interference
            if (distance13 < queryLinkDistance) {
                if (distance23_pair.second < activeChildLinkDistance) {
                    return true;  // return true for interference found
                }
            }

            // dmax lune side2 check
            if ((distance23_pair.second - dmaxCandidateInterferingIndex) >= activeChildLinkDistance) {
                flag_add = false;
            }
        }

        if (flag_add) {
            activeInterferingList.insert(candidateInterferingPivotIndex);
        }
    }

    return false;  // return false for no interference occured by brute force check
}

/**
 * @brief collect potentially interfering fine pivots/exemplar as children of coarse ones
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance link distance from query to activeChildIndex
 * @param activeChildLinkDistance link distance from activeChildIndex to query
 * @param finerLayerIndex
 * @param pivotLayers
 * @param coarseInterferingPivotsList list of interfering coarse pivots
 * @param fineInterferingPivotsList list of interfering fine pivots
 */
void IncrementalBuild_Stage8::collectCandidateFineIndices(unsigned int const queryIndex, unsigned int const activeChildIndex,
                                                          float const queryLinkDistance, float const activeChildLinkDistance,
                                                          int const pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                          tsl::sparse_set<unsigned int> const& coarseInterferingPivotsList,
                                                          tsl::sparse_set<unsigned int>& interferingPivotsList) {
    interferingPivotsList.clear();

    // iterate through each coarse pivot
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarseInterferingPivotsList.begin(); it1 != coarseInterferingPivotsList.end(); it1++) {
        unsigned int const coarseInterferingIndex = (*it1);
        tsl::sparse_set<unsigned int> const& coarseInterferingChildrenList =
            (*pivotLayers)[pivotLayerIndex - 1].get_descendantPivotIndices(coarseInterferingIndex, pivotLayerIndex);

        // consider all children as potentially interfering
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = coarseInterferingChildrenList.begin(); it2 != coarseInterferingChildrenList.end(); it2++) {
            unsigned int const interferingIndex = (*it2);
            bool flag_add = true;

            if (interferingIndex == activeChildIndex) continue;

            // check that all parents of childIndex are in activeCoarseInterferingPivotsList
            tsl::sparse_set<unsigned int> const& fineInterferingParentsList =
                (*pivotLayers)[pivotLayerIndex].get_ancestorPivotIndices(interferingIndex, pivotLayerIndex - 1);
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = fineInterferingParentsList.begin(); it3 != fineInterferingParentsList.end(); it3++) {
                unsigned int const fineInterferingParentIndex = (*it3);
                if (coarseInterferingPivotsList.find(fineInterferingParentIndex) == coarseInterferingPivotsList.end()) {
                    flag_add = false;
                    break;
                }
            }
            if (!flag_add) continue;

            // assisted parent-child check
            flag_add = assistedChildCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, coarseInterferingIndex,
                                          interferingIndex, pivotLayers);

            if (flag_add) {
                interferingPivotsList.insert(interferingIndex);
            }
        }
    }
}

/**
 * @brief check if potentially interfering index is out of range by relationship between link and index parent. Similar to dmax, but for that child in
 * particular
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param coarseInterferingPivot
 * @param fineInterferingPivot
 * @param pivotLayers
 * @return true if fineInterferingPivot can still interfere with the link
 * @return false if out of range of the link
 */
bool IncrementalBuild_Stage8::assistedChildCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                                 float const activeChildLinkDistance, unsigned int const coarseInterferingPivot,
                                                 unsigned int const fineInterferingPivot, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers) {
    bool flag_canInterfere = true;

    // 1: query, 2: activeChild, 3: coarse interfering (parent), 4: fine interfering (child)
    std::pair<bool, float> distance13_pair = pivotLayers->getDistanceIfAvailable(queryIndex, coarseInterferingPivot);
    std::pair<bool, float> distance34_pair = pivotLayers->getDistanceIfAvailable(coarseInterferingPivot, fineInterferingPivot);
    std::pair<bool, float> distance23_pair = pivotLayers->getDistanceIfAvailable(activeChildIndex, coarseInterferingPivot);

    if (distance13_pair.first && distance34_pair.first) {
        if ((distance13_pair.second - distance34_pair.second) >= queryLinkDistance) {
            flag_canInterfere = false;
        }
    }

    if (flag_canInterfere) {
        if (distance23_pair.first && distance34_pair.first) {
            if ((distance23_pair.second - distance34_pair.second) >= activeChildLinkDistance) {
                flag_canInterfere = false;
            }
        }
    }

    return flag_canInterfere;
}

/**
 * @brief brute-force check for interference
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param interferingPivotsList list of potentially interfering pivots
 * @param pivotLayers
 * @return true if interference found
 * @return false if no interference found
 */
bool IncrementalBuild_Stage8::bruteForceInterferenceCheck(unsigned int const queryIndex, unsigned int const activeChildIndex,
                                                          float const queryLinkDistance, float const activeChildLinkDistance,
                                                          tsl::sparse_set<unsigned int> const& interferingPivotsList,
                                                          std::shared_ptr<PivotLayerHierarchy> const& pivotLayers) {
    // test each index for interference
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = interferingPivotsList.begin(); it1 != interferingPivotsList.end(); it1++) {
        unsigned int const interferingPivotIndex = (*it1);

        // brute force interference check
        // could potentially switch this order for less comps...
        float const distance13 = pivotLayers->getDistance(queryIndex, interferingPivotIndex);
        if (distance13 < queryLinkDistance) {
            float const distance23 = pivotLayers->getDistance(activeChildIndex, interferingPivotIndex);
            if (distance23 < activeChildLinkDistance) {
                return true;  // interference found, return true
            }
        }
    }

    return false;  // no interference, becomes a link
}

/**
 *
 *
 *
 *
 *
 *
 *
 *              MULTITHREADING
 *
 *
 *
 *
 *
 *
 */

namespace IncrementalBuild_Stage8_Multithreading {

bool checkLinkForInterference(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex,
                              std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, SparseMatrix& sparseMatrix,
                              SparseMatrixLocal& localSparseMatrix);

void collectCandidateCoarseIndices(unsigned int const queryIndex, int const layerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                   SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix,
                                   tsl::sparse_set<unsigned int> const& activeInterferingPivotsList, RankOrderedList& rankedCandidateList);

bool dmaxCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance, float const activeChildLinkDistance,
               int const layerIndex, int const finerLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, SparseMatrix& sparseMatrix,
               SparseMatrixLocal& localSparseMatrix, RankOrderedList const& candidateInterferingPivots,
               tsl::sparse_set<unsigned int>& activeInterferingList);

void collectCandidateFineIndices(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                 float const activeChildLinkDistance, int const finerLayerIndex,
                                 std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, SparseMatrix& sparseMatrix,
                                 SparseMatrixLocal& localSparseMatrix, tsl::sparse_set<unsigned int> const& activeCoarseInterferingPivotsList,
                                 tsl::sparse_set<unsigned int>& activeFineInterferingPivotsList);

bool assistedChildCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                        float const activeChildLinkDistance, unsigned int const coarseInterferingPivot, unsigned int const fineInterferingPivot,
                        SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix);

bool bruteForceInterferenceCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                 float const activeChildLinkDistance, tsl::sparse_set<unsigned int> const& interferingPivotsList,
                                 SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix);

};  // namespace IncrementalBuild_Stage8_Multithreading

void IncrementalBuild::stage8_multithreading(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int>& finerLayerPivotNeighborsList = queryStruct.finerCoarseNeighborsList;

    // collect activeChildren as children of neighbors
    tsl::sparse_set<unsigned int> activeChildrenList{};
    IncrementalBuild_Stage8::collectNeighborsChildren(queryStruct, _pivotLayers, activeChildrenList);
    std::vector<unsigned int> activeChildrenListVector(activeChildrenList.begin(), activeChildrenList.end());

#pragma omp parallel num_threads(_numThreads)
    {
        // create local sparse matrix and local neighbors list
        SparseMatrixLocal localSparseMatrix = _sparseMatrix->getLocalCopy();
        tsl::sparse_set<unsigned int> localFinerLayerPivotNeighborsList{};

#pragma omp for schedule(static)
        for (int i = 0; i < (int)activeChildrenListVector.size(); i++) {
            unsigned int const activeChildIndex = activeChildrenListVector[i];

            // exhaustive, top-down search for interference
            bool flag_interference = IncrementalBuild_Stage8_Multithreading::checkLinkForInterference(queryStruct, activeChildIndex, _pivotLayers,
                                                                                                      *_sparseMatrix, localSparseMatrix);
            if (!flag_interference) {
                localFinerLayerPivotNeighborsList.insert(activeChildIndex);
            }
        }

#pragma omp critical(coarseNeighborsListUpdate)
        { finerLayerPivotNeighborsList.insert(localFinerLayerPivotNeighborsList.begin(), localFinerLayerPivotNeighborsList.end()); }

#pragma omp critical(sparseMatrixUpdate)
        { _sparseMatrix->updateSparseMatrixWithLocal(localSparseMatrix); }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[8] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[8] += (dEnd_stage - dStart_stage);
    return;
}

/**
 * @brief Exhaustively checks potential link(query,activeChildIndex) for interference in top-down manner
 *
 * @param queryStruct
 * @param activeChildIndex
 * @param finerLayerIndex
 * @param pivotLayers
 * @return true
 * @return false
 */
bool IncrementalBuild_Stage8_Multithreading::checkLinkForInterference(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex,
                                                                      std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                                      SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix) {
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    float const activeChildRadius = (*pivotLayers)[pivotLayerIndex + 1].get_pivotRadius(activeChildIndex);
    float const distance12 = sparseMatrix.getDistance_LocalAndGlobal(queryIndex, activeChildIndex, localSparseMatrix);
    float queryLinkDistance = pivotLayers->_luneRadius(distance12, queryRadius, activeChildRadius);
    float activeChildLinkDistance = pivotLayers->_luneRadius(distance12, activeChildRadius, queryRadius);

    // no lune, return no intererence
    if (queryLinkDistance <= 0 || activeChildLinkDistance <= 0) {
        return false;
    }

    bool flag_interference = false;

    // //------------------------------------------------------------------------------
    // // Work top-down to collect all coarse pivots whose children could interfere
    // //------------------------------------------------------------------------------
    tsl::sparse_set<unsigned int> interferingPivotsList{};
    if (pivotLayerIndex == 0) {
        tsl::sparse_set<unsigned int> const& indices = *(*pivotLayers)[pivotLayerIndex].get_pivotIndices_ptr();
        interferingPivotsList = indices;
    } else {
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>> interferingPivotsMap;
        for (int layerIndex = 0; layerIndex < pivotLayerIndex; layerIndex++) {
            interferingPivotsMap[layerIndex] = tsl::sparse_set<unsigned int>{};

            // created a ranked-list of candidate interfering pivots
            RankOrderedList candidateInterferingPivots;
            if (layerIndex == 0) {
                candidateInterferingPivots = queryStruct.rankOrderedTopLayerPivotsList;
            } else {
                collectCandidateCoarseIndices(queryIndex, layerIndex, pivotLayers, sparseMatrix, localSparseMatrix,
                                              interferingPivotsMap[layerIndex - 1], candidateInterferingPivots);
            }

            // create list of actively interfering pivots from the candidate interfering using dmax
            // return true if any pivots found to be interfering
            flag_interference = dmaxCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, layerIndex, pivotLayerIndex,
                                          pivotLayers, sparseMatrix, localSparseMatrix, candidateInterferingPivots, interferingPivotsMap[layerIndex]);

            if (flag_interference == true) {
                return true;  // interference, this link is no longer to be considered!
            }
        }

        //------------------------------------------------------------------------------
        // Collect potentially interfering fine pivots from children of interfering coarse pivots
        //------------------------------------------------------------------------------
        collectCandidateFineIndices(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, pivotLayerIndex, pivotLayers, sparseMatrix, localSparseMatrix,
                                    interferingPivotsMap[pivotLayerIndex - 1], interferingPivotsList);
    }

    //------------------------------------------------------------------------------
    // Finally, brute force check for interference
    //------------------------------------------------------------------------------
    flag_interference = bruteForceInterferenceCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, interferingPivotsList,
                                                    sparseMatrix, localSparseMatrix);

    return flag_interference;
}

/**
 * @brief Collect potentially interfering fine pivots as children of coarse potentially interfering
 * collectCandidateCoarseIndices(queryIndex, layerIndex, pivotLayers, interferingPivotsMap[layerIndex - 1], candidateInterferingPivots);
 * @param queryIndex
 * @param layerIndex
 * @param pivotLayers
 * @param activeInterferingPivotsList
 * @param rankedCandidateList
 */
void IncrementalBuild_Stage8_Multithreading::collectCandidateCoarseIndices(unsigned int const queryIndex, int const layerIndex,
                                                                           std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                                           SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix,
                                                                           tsl::sparse_set<unsigned int> const& activeInterferingPivotsList,
                                                                           RankOrderedList& rankedCandidateList) {
    rankedCandidateList.clear();

    // iterate through each active pivot in layerIndex-1
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeInterferingPivotsList.begin(); it1 != activeInterferingPivotsList.end(); it1++) {
        unsigned int const activeIndex = (*it1);
        tsl::sparse_set<unsigned int> const& activeIndexChildrenList =
            (*pivotLayers)[layerIndex - 1].get_descendantPivotIndices(activeIndex, layerIndex);

        // consider all children of activeIndex as potential candidate
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activeIndexChildrenList.begin(); it2 != activeIndexChildrenList.end(); it2++) {
            unsigned int const childIndex = (*it2);
            bool flag_add = true;

            // check that all parents of childIndex are in activeInterferingPivotsList
            tsl::sparse_set<unsigned int> const& childIndexParentsList =
                (*pivotLayers)[layerIndex].get_ancestorPivotIndices(childIndex, layerIndex - 1);
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = childIndexParentsList.begin(); it3 != childIndexParentsList.end(); it3++) {
                unsigned int const childParentIndex = (*it3);

                if (activeInterferingPivotsList.find(childParentIndex) == activeInterferingPivotsList.end()) {
                    flag_add = false;
                    break;
                }
            }

            // child becomes potentially interfering pivot
            if (flag_add) {
                float const distance = sparseMatrix.getDistance_LocalAndGlobal(queryIndex, childIndex, localSparseMatrix);
                rankedCandidateList.add(childIndex, distance);
            }
        }
    }
    rankedCandidateList.sort();
}

/**
 * @brief Given list of potentially interfering pivots, use dmax to tell if any children within their domain may interfere with link.
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param layerIndex
 * @param finerLayerIndex
 * @param pivotLayers
 * @param candidateInterferingPivots
 * @param activeInterferingList
 * @return true
 * @return false
 */
bool IncrementalBuild_Stage8_Multithreading::dmaxCheck(unsigned int const queryIndex, unsigned int const activeChildIndex,
                                                       float const queryLinkDistance, float const activeChildLinkDistance, int const layerIndex,
                                                       int const pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                       SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix,
                                                       RankOrderedList const& candidateInterferingPivots,
                                                       tsl::sparse_set<unsigned int>& activeInterferingList) {
    activeInterferingList.clear();

    // create list of actively interfering pivots
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = candidateInterferingPivots.list.begin(); it1 != candidateInterferingPivots.list.end(); it1++) {
        unsigned int const candidateInterferingPivotIndex = (*it1).index;
        float const distance13 = (*it1).distance;
        bool flag_add = true;

        // // check lune side1 against maxDmax, the max of all possible parent-child distances
        float const dmaxMaxLayer = (*pivotLayers)[layerIndex].get_dmaxMaximumOfMaxChildDistance(pivotLayerIndex);
        if ((distance13 - dmaxMaxLayer) >= queryLinkDistance) {
            break;  // break since ranked list
        }

        // check lune side1 against dmax
        float const dmaxCandidateInterferingIndex =
            (*pivotLayers)[layerIndex].get_dmaxMaxChildDistance(candidateInterferingPivotIndex, pivotLayerIndex);
        if ((distance13 - dmaxCandidateInterferingIndex) >= queryLinkDistance) {
            flag_add = false;
            continue;  // no children of this one can interfere, next one!
        }

        // if distance available, lune side 2 with dmax and brute force check
        std::pair<bool, float> distance23_pair =
            sparseMatrix.getDistanceIfAvailable_LocalAndGlobal(activeChildIndex, candidateInterferingPivotIndex, localSparseMatrix);
        if (distance23_pair.first) {
            // brute force check for pivot interference
            if (distance13 < queryLinkDistance) {
                if (distance23_pair.second < activeChildLinkDistance) {
                    return true;  // return true for interference found
                }
            }

            // dmax lune side2 check
            if ((distance23_pair.second - dmaxCandidateInterferingIndex) >= activeChildLinkDistance) {
                flag_add = false;
            }
        }

        if (flag_add) {
            activeInterferingList.insert(candidateInterferingPivotIndex);
        }
    }

    return false;  // return false for no interference occured by brute force check
}

/**
 * @brief collect potentially interfering fine pivots/exemplar as children of coarse ones
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance link distance from query to activeChildIndex
 * @param activeChildLinkDistance link distance from activeChildIndex to query
 * @param finerLayerIndex
 * @param pivotLayers
 * @param coarseInterferingPivotsList list of interfering coarse pivots
 * @param fineInterferingPivotsList list of interfering fine pivots
 */
void IncrementalBuild_Stage8_Multithreading::collectCandidateFineIndices(unsigned int const queryIndex, unsigned int const activeChildIndex,
                                                                         float const queryLinkDistance, float const activeChildLinkDistance,
                                                                         int const pivotLayerIndex,
                                                                         std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                                         SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix,
                                                                         tsl::sparse_set<unsigned int> const& coarseInterferingPivotsList,
                                                                         tsl::sparse_set<unsigned int>& interferingPivotsList) {
    interferingPivotsList.clear();

    // iterate through each coarse pivot
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarseInterferingPivotsList.begin(); it1 != coarseInterferingPivotsList.end(); it1++) {
        unsigned int const coarseInterferingIndex = (*it1);
        tsl::sparse_set<unsigned int> const& coarseInterferingChildrenList =
            (*pivotLayers)[pivotLayerIndex - 1].get_descendantPivotIndices(coarseInterferingIndex, pivotLayerIndex);

        // consider all children as potentially interfering
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = coarseInterferingChildrenList.begin(); it2 != coarseInterferingChildrenList.end(); it2++) {
            unsigned int const interferingIndex = (*it2);
            bool flag_add = true;

            if (interferingIndex == activeChildIndex) continue;

            // check that all parents of childIndex are in activeCoarseInterferingPivotsList
            tsl::sparse_set<unsigned int> const& fineInterferingParentsList =
                (*pivotLayers)[pivotLayerIndex].get_ancestorPivotIndices(interferingIndex, pivotLayerIndex - 1);
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = fineInterferingParentsList.begin(); it3 != fineInterferingParentsList.end(); it3++) {
                unsigned int const fineInterferingParentIndex = (*it3);
                if (coarseInterferingPivotsList.find(fineInterferingParentIndex) == coarseInterferingPivotsList.end()) {
                    flag_add = false;
                    break;
                }
            }
            if (!flag_add) continue;

            // assisted parent-child check
            flag_add = assistedChildCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, coarseInterferingIndex,
                                          interferingIndex, sparseMatrix, localSparseMatrix);

            if (flag_add) {
                interferingPivotsList.insert(interferingIndex);
            }
        }
    }
}

/**
 * @brief check if potentially interfering index is out of range by relationship between link and index parent. Similar to dmax, but for that child in
 * particular
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param coarseInterferingPivot
 * @param fineInterferingPivot
 * @param pivotLayers
 * @return true if fineInterferingPivot can still interfere with the link
 * @return false if out of range of the link
 */
bool IncrementalBuild_Stage8_Multithreading::assistedChildCheck(unsigned int const queryIndex, unsigned int const activeChildIndex,
                                                                float const queryLinkDistance, float const activeChildLinkDistance,
                                                                unsigned int const coarseInterferingPivot, unsigned int const fineInterferingPivot,
                                                                SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix) {
    bool flag_canInterfere = true;

    // 1: query, 2: activeChild, 3: coarse interfering (parent), 4: fine interfering (child)
    std::pair<bool, float> distance13_pair =
        sparseMatrix.getDistanceIfAvailable_LocalAndGlobal(queryIndex, coarseInterferingPivot, localSparseMatrix);
    std::pair<bool, float> distance34_pair =
        sparseMatrix.getDistanceIfAvailable_LocalAndGlobal(coarseInterferingPivot, fineInterferingPivot, localSparseMatrix);
    std::pair<bool, float> distance23_pair =
        sparseMatrix.getDistanceIfAvailable_LocalAndGlobal(activeChildIndex, coarseInterferingPivot, localSparseMatrix);

    if (distance13_pair.first && distance34_pair.first) {
        if ((distance13_pair.second - distance34_pair.second) >= queryLinkDistance) {
            flag_canInterfere = false;
        }
    }

    if (flag_canInterfere) {
        if (distance23_pair.first && distance34_pair.first) {
            if ((distance23_pair.second - distance34_pair.second) >= activeChildLinkDistance) {
                flag_canInterfere = false;
            }
        }
    }

    return flag_canInterfere;
}

/**
 * @brief brute-force check for interference
 *
 * @param queryIndex
 * @param activeChildIndex
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param interferingPivotsList list of potentially interfering pivots
 * @param pivotLayers
 * @return true if interference found
 * @return false if no interference found
 */
bool IncrementalBuild_Stage8_Multithreading::bruteForceInterferenceCheck(unsigned int const queryIndex, unsigned int const activeChildIndex,
                                                                         float const queryLinkDistance, float const activeChildLinkDistance,
                                                                         tsl::sparse_set<unsigned int> const& interferingPivotsList,
                                                                         SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix) {
    // test each index for interference
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = interferingPivotsList.begin(); it1 != interferingPivotsList.end(); it1++) {
        unsigned int const interferingPivotIndex = (*it1);

        // brute force interference check
        // could potentially switch this order for less comps...
        float const distance13 = sparseMatrix.getDistance_LocalAndGlobal(queryIndex, interferingPivotIndex, localSparseMatrix);
        if (distance13 < queryLinkDistance) {
            float const distance23 = sparseMatrix.getDistance_LocalAndGlobal(activeChildIndex, interferingPivotIndex, localSparseMatrix);
            if (distance23 < activeChildLinkDistance) {
                return true;  // interference found, return true
            }
        }
    }

    return false;  // no interference, becomes a link
}

/**
 *
 *
 *
 *
 *
 *
 *
 *              BRUTE FORCE
 *
 *
 *
 *
 *
 *
 *
 *
 */

/**
 * @brief Find fine pivots that are virtual neighbors of Q
 *
 * @param queryStruct
 */
void IncrementalBuild::stage8_bruteForce(QueryStructIncrementalBuild& queryStruct) {
    // GET ALL STRUCTURES FROM QIB
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    int const finerPivotLayerIndex = pivotLayerIndex + 1;
    tsl::sparse_set<unsigned int>& finerCoarseNeighborsList = queryStruct.finerCoarseNeighborsList;

    // collect children of all neighbors as potential virtual neighbors
    tsl::sparse_set<unsigned int> candidateNeighborList{};
    // IncrementalBuild_Stage8::collectNeighborsChildren(queryStruct, _pivotLayers, candidateNeighborList);
    tsl::sparse_set<unsigned int> const& finePivotIndices = *(*_pivotLayers)[finerPivotLayerIndex].get_pivotIndices_ptr();
    candidateNeighborList = finePivotIndices;

    // perform brute force neighbors validation. check all pivots for interference
    tsl::sparse_set<unsigned int> const& pivotIndices = *(*_pivotLayers)[pivotLayerIndex].get_pivotIndices_ptr();
    tsl::sparse_set<unsigned int>::const_iterator it2, it3;
    for (it2 = candidateNeighborList.begin(); it2 != candidateNeighborList.end(); it2++) {
        unsigned int const finePivotIndex = (*it2);
        float const finePivotRadius = (*_pivotLayers)[finerPivotLayerIndex].get_pivotRadius(finePivotIndex);
        float const distance12 = _sparseMatrix->getDistance(queryIndex, finePivotIndex);
        bool flag_neighbor = true;

        // if no lune, automatic neighbor
        if (distance12 - (2 * queryRadius + finePivotRadius) <= 0 || distance12 - (queryRadius + 2 * finePivotRadius) <= 0) {
            finerCoarseNeighborsList.insert(finePivotIndex);
            continue;
        }

        // check for interference by ALL indices
        for (it3 = pivotIndices.begin(); it3 != pivotIndices.end(); it3++) {
            unsigned int const interferingPivotIndex = (*it3);
            float const distance13 = _sparseMatrix->getDistance(queryIndex, interferingPivotIndex);
            if (distance13 < _luneRadius(distance12, queryRadius, finePivotRadius)) {
                float const distance23 = _sparseMatrix->getDistance(finePivotIndex, interferingPivotIndex);
                if (distance23 < _luneRadius(distance12, finePivotRadius, queryRadius)) {
                    flag_neighbor = false;
                    break;
                }
            }
        }

        if (flag_neighbor) {
            finerCoarseNeighborsList.insert(finePivotIndex);
        }
    }

    return;
}
