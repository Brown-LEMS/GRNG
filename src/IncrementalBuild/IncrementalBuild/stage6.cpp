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
// 2022-05-31
#include "incremental-build.hpp"

//  Helper Functions for OfflineBuild::stage6_multithreading()

bool IB0_checkLinkForInterference_multithreading(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex, int const finerLayerIndex,
                                             std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, SparseMatrix& sparseMatrix,
                                             SparseMatrixLocal& localSparseMatrix);
void IB1_collectCandidateCoarseIndices_multithreading(unsigned int const queryIndex, int const layerIndex,
                                                     std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                     tsl::sparse_set<unsigned int> const& activeInterferingPivotsList,
                                                     RankOrderedList& rankedCandidateList, SparseMatrix& sparseMatrix,
                                                     SparseMatrixLocal& localSparseMatrix);
bool IB2_dmaxCheck_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                 float const activeChildLinkDistance, int const layerIndex, int const finerLayerIndex,
                                 std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, RankOrderedList const& candidateInterferingPivots,
                                 tsl::sparse_set<unsigned int>& activeInterferingList, SparseMatrix& sparseMatrix,
                                 SparseMatrixLocal& localSparseMatrix);
void IB3_collectCandidateFineIndices_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                                   float const activeChildLinkDistance, int const finerLayerIndex,
                                                   std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                   tsl::sparse_set<unsigned int> const& coarseInterferingPivotsList,
                                                   tsl::sparse_set<unsigned int>& fineInterferingPivotsList, SparseMatrix& sparseMatrix,
                                                   SparseMatrixLocal& localSparseMatrix);
bool IB3a_assistedChildCheck_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                           float const activeChildLinkDistance, unsigned int const coarseInterferingPivot,
                                           unsigned int const fineInterferingPivot, SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix);
bool IB4_bruteForceInterferenceCheck_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                                   float const activeChildLinkDistance, tsl::sparse_set<unsigned int> const& interferingPivotsList,
                                                   SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix);

/**
 * @brief RNG Link Verification: exhaustive check for interfence using multithreading
 *
 * @param queryStruct
 */
void IncrementalBuild::stage6_multithreading(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int> const& activeChildrenList = queryStruct.activeChildrenList;
    tsl::sparse_set<unsigned int>& finerLayerPivotNeighborsList = queryStruct.pivotLayerNeighborsListMap[finerLayerIndex];

    // vector for #pragma loop
    std::vector<unsigned int> activeChildrenListVector(activeChildrenList.begin(), activeChildrenList.end());

#pragma omp parallel num_threads(_numThreads)
    {
        // create local sparse matrix and local neighbors list
        SparseMatrixLocal localSparseMatrix = _sparseMatrix->getLocalCopy();
        tsl::sparse_set<unsigned int> localFinerPivotLayerNeighborsList;

#pragma omp for schedule(static)
        for (int i = 0; i < (int)activeChildrenListVector.size(); i++) {
            unsigned int const activeChildIndex = activeChildrenListVector[i];

            // exhaustive, top-down search for interference
            bool flag_interference = IB0_checkLinkForInterference_multithreading(queryStruct, activeChildIndex, finerLayerIndex, _pivotLayers,
                                                                             *_sparseMatrix, localSparseMatrix);
            if (!flag_interference) {
                localFinerPivotLayerNeighborsList.insert(activeChildIndex);
            }
        }

#pragma omp critical(pivotLayerNeighborsListUpdate)
        { finerLayerPivotNeighborsList.insert(localFinerPivotLayerNeighborsList.begin(), localFinerPivotLayerNeighborsList.end()); }

#pragma omp critical(sparseMatrixUpdate)
        { _sparseMatrix->updateSparseMatrixWithLocal(localSparseMatrix); }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[6] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[6] += (dEnd_stage - dStart_stage);

    return;
}

/**
 * @brief
 *
 * @param queryStruct
 * @param activeChildIndex
 * @param finerLayerIndex
 * @param pivotLayers
 * @param sparseMatrix
 * @param localSparseMatrix
 * @return true
 * @return false
 */
bool IB0_checkLinkForInterference_multithreading(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex, int const finerLayerIndex,
                                             std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, SparseMatrix& sparseMatrix,
                                             SparseMatrixLocal& localSparseMatrix) {
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    float const activeChildRadius = (*pivotLayers)[finerLayerIndex].get_pivotRadius(activeChildIndex);
    float const distance12 = sparseMatrix.getDistance_LocalAndGlobal(queryIndex, activeChildIndex, localSparseMatrix);
    float queryLinkDistance = pivotLayers->_luneRadius(distance12, queryRadius, activeChildRadius);
    float activeChildLinkDistance = pivotLayers->_luneRadius(distance12, activeChildRadius, queryRadius);
    bool flag_interference = false;

    //------------------------------------------------------------------------------
    // Work top-down to collect all coarse pivots whose children could interfere
    //------------------------------------------------------------------------------
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> interferingPivotsMap;
    for (int layerIndex = 0; layerIndex < finerLayerIndex; layerIndex++) {
        interferingPivotsMap[layerIndex] = tsl::sparse_set<unsigned int>{};

        // created a ranked-list of candidate interfering pivots
        RankOrderedList candidateInterferingPivots;
        if (layerIndex == 0) {
            candidateInterferingPivots = queryStruct.rankOrderedTopLayerPivotsList;
        } else {
            IB1_collectCandidateCoarseIndices_multithreading(queryIndex, layerIndex, pivotLayers, interferingPivotsMap[layerIndex - 1],
                                                            candidateInterferingPivots, sparseMatrix, localSparseMatrix);
        }

        // create list of actively interfering pivots from the candidate interfering using dmax
        // return true if any pivots found to be interfering
        flag_interference =
            IB2_dmaxCheck_multithreading(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, layerIndex, finerLayerIndex,
                                        pivotLayers, candidateInterferingPivots, interferingPivotsMap[layerIndex], sparseMatrix, localSparseMatrix);

        if (flag_interference) return true;  // interference, this link is no longer to be considered!
    }

    //------------------------------------------------------------------------------
    // Collect potentially interfering fine pivots from children of interfering coarse pivots
    //------------------------------------------------------------------------------
    tsl::sparse_set<unsigned int> fineInterferingPivotsList;
    IB3_collectCandidateFineIndices_multithreading(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, finerLayerIndex,
                                                  pivotLayers, interferingPivotsMap[finerLayerIndex - 1], fineInterferingPivotsList, sparseMatrix,
                                                  localSparseMatrix);

    //------------------------------------------------------------------------------
    // Finally, brute force check for interference
    //------------------------------------------------------------------------------
    flag_interference = IB4_bruteForceInterferenceCheck_multithreading(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance,
                                                                      fineInterferingPivotsList, sparseMatrix, localSparseMatrix);

    return flag_interference;
}

/**
 * @brief Collect potentially interfering fine pivots as children of coarse potentially interfering
 *
 * @param queryIndex
 * @param layerIndex
 * @param pivotLayers
 * @param activeInterferingPivotsList
 * @param rankedCandidateList
 * @param sparseMatrix
 * @param localSparseMatrix
 */
void IB1_collectCandidateCoarseIndices_multithreading(unsigned int const queryIndex, int const layerIndex,
                                                     std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                     tsl::sparse_set<unsigned int> const& activeInterferingPivotsList,
                                                     RankOrderedList& rankedCandidateList, SparseMatrix& sparseMatrix,
                                                     SparseMatrixLocal& localSparseMatrix) {
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
    return;
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
 * @param sparseMatrix
 * @param localSparseMatrix
 * @return true
 * @return false
 */
bool IB2_dmaxCheck_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                 float const activeChildLinkDistance, int const layerIndex, int const finerLayerIndex,
                                 std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, RankOrderedList const& candidateInterferingPivots,
                                 tsl::sparse_set<unsigned int>& activeInterferingList, SparseMatrix& sparseMatrix,
                                 SparseMatrixLocal& localSparseMatrix) {
    activeInterferingList.clear();

    // create list of actively interfering pivots
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = candidateInterferingPivots.list.begin(); it1 != candidateInterferingPivots.list.end(); it1++) {
        unsigned int const candidateInterferingPivotIndex = (*it1).index;
        float const distance13 = (*it1).distance;
        bool flag_add = true;

        // check lune side1 against maxDmax, the max of all possible parent-child distances
        float const dmaxMaxLayer = (*pivotLayers)[layerIndex].get_dmaxMaximumOfMaxChildDistance(finerLayerIndex);
        if ((distance13 - dmaxMaxLayer) >= queryLinkDistance) {
            break;  // break since ranked list
        }

        // check lune side1 against dmax
        float const dmaxCandidateInterferingIndex =
            (*pivotLayers)[layerIndex].get_dmaxMaxChildDistance(candidateInterferingPivotIndex, finerLayerIndex);
        if ((distance13 - dmaxCandidateInterferingIndex) >= queryLinkDistance) {
            flag_add = false;
            continue;  // no children of this one can interfere, next one!
        }

        // if distance available, lune side 2 with dmax and brute force check
        std::pair<bool, float> const distance23_pair =
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
 * @param queryLinkDistance
 * @param activeChildLinkDistance
 * @param finerLayerIndex
 * @param pivotLayers
 * @param activeCoarseInterferingPivotsList
 * @param activeFineInterferingPivotsList
 * @param sparseMatrix
 * @param localSparseMatrix
 */
void IB3_collectCandidateFineIndices_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                                   float const activeChildLinkDistance, int const finerLayerIndex,
                                                   std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                   tsl::sparse_set<unsigned int> const& coarseInterferingPivotsList,
                                                   tsl::sparse_set<unsigned int>& fineInterferingPivotsList, SparseMatrix& sparseMatrix,
                                                   SparseMatrixLocal& localSparseMatrix) {
    fineInterferingPivotsList.clear();

    // iterate through each coarse pivot
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarseInterferingPivotsList.begin(); it1 != coarseInterferingPivotsList.end(); it1++) {
        unsigned int const coarseInterferingIndex = (*it1);
        tsl::sparse_set<unsigned int> const& coarseInterferingChildrenList =
            (*pivotLayers)[finerLayerIndex - 1].get_descendantPivotIndices(coarseInterferingIndex, finerLayerIndex);

        // consider all children as potentially interfering
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = coarseInterferingChildrenList.begin(); it2 != coarseInterferingChildrenList.end(); it2++) {
            unsigned int const fineInterferingIndex = (*it2);
            bool flag_add = true;

            if (fineInterferingIndex == activeChildIndex) continue;

            // check that all parents of childIndex are in activeCoarseInterferingPivotsList
            tsl::sparse_set<unsigned int> const& fineInterferingParentsList =
                (*pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(fineInterferingIndex, finerLayerIndex - 1);
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
            flag_add = IB3a_assistedChildCheck_multithreading(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance,
                                                             coarseInterferingIndex, fineInterferingIndex, sparseMatrix, localSparseMatrix);

            if (flag_add) {
                fineInterferingPivotsList.insert(fineInterferingIndex);
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
 * @param sparseMatrix
 * @param localSparseMatrix
 * @return true
 * @return false
 */
bool IB3a_assistedChildCheck_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                           float const activeChildLinkDistance, unsigned int const coarseInterferingPivot,
                                           unsigned int const fineInterferingPivot, SparseMatrix& sparseMatrix,
                                           SparseMatrixLocal& localSparseMatrix) {
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
 * @param interferingPivotsList
 * @param sparseMatrix
 * @param localSparseMatrix
 * @return true
 * @return false
 */
bool IB4_bruteForceInterferenceCheck_multithreading(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                                   float const activeChildLinkDistance, tsl::sparse_set<unsigned int> const& interferingPivotsList,
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

// ===================================================================================================================================================
//
//
//
//
//
//
//
//
//
//
//
//                                      MULTITHREAD VS NORMAL SEPARATION
//
//
//
//
//
//
//
//
//
//
//
// ===================================================================================================================================================

bool IB0_checkLinkForInterference(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex, int const finerLayerIndex,
                              std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

void IB1_collectCandidateCoarseIndices(unsigned int const queryIndex, int const layerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                      tsl::sparse_set<unsigned int> const& activeInterferingPivotsList, RankOrderedList& rankedCandidateList);

bool IB2_dmaxCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                  float const activeChildLinkDistance, int const layerIndex, int const finerLayerIndex,
                  std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, RankOrderedList const& candidateInterferingPivots,
                  tsl::sparse_set<unsigned int>& activeInterferingList);

void IB3_collectCandidateFineIndices(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                    float const activeChildLinkDistance, int const finerLayerIndex,
                                    std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                    tsl::sparse_set<unsigned int> const& activeCoarseInterferingPivotsList,
                                    tsl::sparse_set<unsigned int>& activeFineInterferingPivotsList);

bool IB3a_assistedChildCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                            float const activeChildLinkDistance, unsigned int const coarseInterferingPivot, unsigned int const fineInterferingPivot,
                            std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

bool IB4_bruteForceInterferenceCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                    float const activeChildLinkDistance, tsl::sparse_set<unsigned int> const& interferingPivotsList,
                                    std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

/**
 * @brief RNG Link Verification: exhaustive check for interfence
 *
 * @param queryStruct
 */
void IncrementalBuild::stage6(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int> const& activeChildrenList = queryStruct.activeChildrenList;
    tsl::sparse_set<unsigned int>& finerLayerPivotNeighborsList = queryStruct.pivotLayerNeighborsListMap[finerLayerIndex];

    // top-down, exhaustive check for interference of link(query,activeChild)
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeChildrenList.begin(); it1 != activeChildrenList.end(); it1++) {
        unsigned int const activeChildIndex = (*it1);

        // exhaustive, top-down search for interference
        bool flag_interference = IB0_checkLinkForInterference(queryStruct, activeChildIndex, finerLayerIndex, _pivotLayers);
        if (!flag_interference) {
            finerLayerPivotNeighborsList.insert(activeChildIndex);
        }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[6] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[6] += (dEnd_stage - dStart_stage);
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
bool IB0_checkLinkForInterference(QueryStructIncrementalBuild& queryStruct, unsigned int const activeChildIndex, int const finerLayerIndex,
                              std::shared_ptr<PivotLayerHierarchy> const& pivotLayers) {
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    float const activeChildRadius = (*pivotLayers)[finerLayerIndex].get_pivotRadius(activeChildIndex);
    float const distance12 = pivotLayers->getDistance(queryIndex, activeChildIndex);
    float queryLinkDistance = pivotLayers->_luneRadius(distance12, queryRadius, activeChildRadius);
    float activeChildLinkDistance = pivotLayers->_luneRadius(distance12, activeChildRadius, queryRadius);
    bool flag_interference = false;

    //------------------------------------------------------------------------------
    // Work top-down to collect all coarse pivots whose children could interfere
    //------------------------------------------------------------------------------
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> interferingPivotsMap;
    for (int layerIndex = 0; layerIndex < finerLayerIndex; layerIndex++) {
        interferingPivotsMap[layerIndex] = tsl::sparse_set<unsigned int>{};

        // created a ranked-list of candidate interfering pivots
        RankOrderedList candidateInterferingPivots;
        if (layerIndex == 0) {
            candidateInterferingPivots = queryStruct.rankOrderedTopLayerPivotsList;
        } else {
            IB1_collectCandidateCoarseIndices(queryIndex, layerIndex, pivotLayers, interferingPivotsMap[layerIndex - 1], candidateInterferingPivots);
        }

        // create list of actively interfering pivots from the candidate interfering using dmax
        // return true if any pivots found to be interfering
        flag_interference = IB2_dmaxCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, layerIndex, finerLayerIndex,
                                         pivotLayers, candidateInterferingPivots, interferingPivotsMap[layerIndex]);

        if (flag_interference) return true;  // interference, this link is no longer to be considered!
    }

    //------------------------------------------------------------------------------
    // Collect potentially interfering fine pivots from children of interfering coarse pivots
    //------------------------------------------------------------------------------
    tsl::sparse_set<unsigned int> fineInterferingPivotsList;
    IB3_collectCandidateFineIndices(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, finerLayerIndex, pivotLayers,
                                   interferingPivotsMap[finerLayerIndex - 1], fineInterferingPivotsList);

    //------------------------------------------------------------------------------
    // Finally, brute force check for interference
    //------------------------------------------------------------------------------
    flag_interference = IB4_bruteForceInterferenceCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance,
                                                       fineInterferingPivotsList, pivotLayers);

    return flag_interference;
}

/**
 * @brief Collect potentially interfering fine pivots as children of coarse potentially interfering
 *
 * @param queryIndex
 * @param layerIndex
 * @param pivotLayers
 * @param activeInterferingPivotsList
 * @param rankedCandidateList
 */
void IB1_collectCandidateCoarseIndices(unsigned int const queryIndex, int const layerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                      tsl::sparse_set<unsigned int> const& activeInterferingPivotsList, RankOrderedList& rankedCandidateList) {
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
bool IB2_dmaxCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                  float const activeChildLinkDistance, int const layerIndex, int const finerLayerIndex,
                  std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, RankOrderedList const& candidateInterferingPivots,
                  tsl::sparse_set<unsigned int>& activeInterferingList) {
    activeInterferingList.clear();

    // create list of actively interfering pivots
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = candidateInterferingPivots.list.begin(); it1 != candidateInterferingPivots.list.end(); it1++) {
        unsigned int const candidateInterferingPivotIndex = (*it1).index;
        float const distance13 = (*it1).distance;
        bool flag_add = true;

        // check lune side1 against maxDmax, the max of all possible parent-child distances
        float const dmaxMaxLayer = (*pivotLayers)[layerIndex].get_dmaxMaximumOfMaxChildDistance(finerLayerIndex);
        if ((distance13 - dmaxMaxLayer) >= queryLinkDistance) {
            break;  // break since ranked list
        }

        // check lune side1 against dmax
        float const dmaxCandidateInterferingIndex =
            (*pivotLayers)[layerIndex].get_dmaxMaxChildDistance(candidateInterferingPivotIndex, finerLayerIndex);
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
void IB3_collectCandidateFineIndices(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                    float const activeChildLinkDistance, int const finerLayerIndex,
                                    std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                    tsl::sparse_set<unsigned int> const& coarseInterferingPivotsList,
                                    tsl::sparse_set<unsigned int>& fineInterferingPivotsList) {
    fineInterferingPivotsList.clear();

    // iterate through each coarse pivot
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarseInterferingPivotsList.begin(); it1 != coarseInterferingPivotsList.end(); it1++) {
        unsigned int const coarseInterferingIndex = (*it1);
        tsl::sparse_set<unsigned int> const& coarseInterferingChildrenList =
            (*pivotLayers)[finerLayerIndex - 1].get_descendantPivotIndices(coarseInterferingIndex, finerLayerIndex);

        // consider all children as potentially interfering
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = coarseInterferingChildrenList.begin(); it2 != coarseInterferingChildrenList.end(); it2++) {
            unsigned int const fineInterferingIndex = (*it2);
            bool flag_add = true;

            if (fineInterferingIndex == activeChildIndex) continue;

            // check that all parents of childIndex are in activeCoarseInterferingPivotsList
            tsl::sparse_set<unsigned int> const& fineInterferingParentsList =
                (*pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(fineInterferingIndex, finerLayerIndex - 1);
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
            flag_add = IB3a_assistedChildCheck(queryIndex, activeChildIndex, queryLinkDistance, activeChildLinkDistance, coarseInterferingIndex,
                                              fineInterferingIndex, pivotLayers);

            if (flag_add) {
                fineInterferingPivotsList.insert(fineInterferingIndex);
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
bool IB3a_assistedChildCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                            float const activeChildLinkDistance, unsigned int const coarseInterferingPivot, unsigned int const fineInterferingPivot,
                            std::shared_ptr<PivotLayerHierarchy> const& pivotLayers) {
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
bool IB4_bruteForceInterferenceCheck(unsigned int const queryIndex, unsigned int const activeChildIndex, float const queryLinkDistance,
                                    float const activeChildLinkDistance, tsl::sparse_set<unsigned int> const& interferingPivotsList,
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
 * @brief Stage 6 by brute force using max number of threads available
 *
 * @param queryStruct
 */
void IncrementalBuild::stage6_bruteForce(QueryStructIncrementalBuild& queryStruct) {
    // GET ALL STRUCTURES FROM QIB
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserPivotLayerIndex = queryStruct.coarserLayerIndex;
    int const finerPivotLayerIndex = coarserPivotLayerIndex + 1;

    tsl::sparse_set<unsigned int> const& activeChildrenList = queryStruct.activeChildrenList;
    tsl::sparse_set<unsigned int>& finerLayerPivotNeighborList = queryStruct.pivotLayerNeighborsListMap[finerPivotLayerIndex];

    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeChildrenList.begin(); it1 != activeChildrenList.end(); it1++) {
        unsigned int const activeChildIndex = (*it1);
        float const activeChildRadius = (*_pivotLayers)[finerPivotLayerIndex].get_pivotRadius(activeChildIndex);
        float const distance12 = _sparseMatrix->getDistance(queryIndex, activeChildIndex);
        bool flag_isNeighbor = true;

        // IF NO LUNE, CANNOT BE INTERFERED WITH, BECOMES A PIVOT
        if (distance12 - (2 * queryRadius + activeChildRadius) <= 0 || distance12 - (queryRadius + 2 * activeChildRadius) <= 0) {
            finerLayerPivotNeighborList.insert(activeChildIndex);
            continue;
        }

        // CHECK ALL PIVOTS AS POTENTIALLY INTERFERING
        tsl::sparse_set<unsigned int> const& pivotInterferenceIndices = *((*_pivotLayers)[finerPivotLayerIndex].get_pivotIndices_ptr());
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = pivotInterferenceIndices.begin(); it2 != pivotInterferenceIndices.end(); it2++) {
            unsigned int const interferingPivotIndex = (*it2);
            if (interferingPivotIndex == queryIndex || interferingPivotIndex == activeChildIndex) continue;

            // BRUTE FORCE INTERFERENCE CHECK
            float const distance13 = _sparseMatrix->getDistance(queryIndex, interferingPivotIndex);
            if (distance13 < distance12 - (2 * queryRadius + activeChildRadius)) {
                float const distance23 = _sparseMatrix->getDistance(activeChildIndex, interferingPivotIndex);
                if (distance23 < distance12 - (queryRadius + 2 * activeChildRadius)) {
                    flag_isNeighbor = false;
                    break;
                }
            }
        }

        // IF NOT INTERFERED WITH, ADD AS NEIGHBOR
        if (flag_isNeighbor) {
            finerLayerPivotNeighborList.insert(activeChildIndex);
        }
    }

    // sort neighbors to always have them in same order despite threads
    //std::vector<unsigned int> finerPivotLayerNeighborsCopy(finerLayerPivotNeighborList->begin(), finerLayerPivotNeighborList->end());
    //std::sort(finerPivotLayerNeighborsCopy.begin(), finerPivotLayerNeighborsCopy.end());
    //finerLayerPivotNeighborList->clear();
    //finerLayerPivotNeighborList->insert(finerPivotLayerNeighborsCopy.begin(), finerPivotLayerNeighborsCopy.end());

    return;
}

// /**
//  * @brief Stage 6 by brute force using max number of threads available
//  *
//  * @param queryStruct
//  */
// void IncrementalBuild::stage6_bruteForce(QueryStructIncrementalBuild& queryStruct) {
//     // GET ALL STRUCTURES FROM QIB
//     unsigned int const queryIndex = queryStruct.queryIndex;
//     float const queryRadius = queryStruct.queryRadius;
//     int const coarserPivotLayerIndex = queryStruct.coarserLayerIndex;
//     int const finerPivotLayerIndex = coarserPivotLayerIndex + 1;

//     tsl::sparse_set<unsigned int> const& activeChildrenList = queryStruct.activeChildrenList;
//     tsl::sparse_set<unsigned int>* finerLayerPivotNeighborList = &queryStruct.pivotLayerNeighborsListMap[finerPivotLayerIndex];
//     std::vector<unsigned int> activeChildrenListVector(activeChildrenList.begin(), activeChildrenList.end());

// #pragma omp parallel num_threads(_numThreads)
//     {
//         // create local sparse matrix and local neighbors list
//         SparseMatrixLocal localSparseMatrix;
//         tsl::sparse_set<unsigned int> localFinerPivotLayerNeighborsList;

// #pragma omp for schedule(static)
//         for (int i = 0; i < (int)activeChildrenListVector.size(); i++) {
//             // potential neighbor information
//             unsigned int const activeChildIndex = activeChildrenListVector[i];
//             float const activeChildRadius = (*_pivotLayers)[finerPivotLayerIndex].get_pivotRadius(activeChildIndex);
//             float const distance12 = _sparseMatrix->getDistance_LocalAndGlobal(queryIndex, activeChildIndex, localSparseMatrix);
//             bool flag_isNeighbor = true;

//             // IF NO LUNE, CANNOT BE INTERFERED WITH, BECOMES A PIVOT
//             if (distance12 - (2 * queryRadius + activeChildRadius) <= 0 || distance12 - (queryRadius + 2 * activeChildRadius) <= 0) {
//                 localFinerPivotLayerNeighborsList.insert(activeChildIndex);
//                 continue;
//             }

//             // CHECK ALL PIVOTS AS POTENTIALLY INTERFERING
//             tsl::sparse_set<unsigned int> const& pivotInterferenceIndices = *((*_pivotLayers)[finerPivotLayerIndex].get_pivotIndices_ptr());
//             tsl::sparse_set<unsigned int>::const_iterator it2;
//             for (it2 = pivotInterferenceIndices.begin(); it2 != pivotInterferenceIndices.end(); it2++) {
//                 unsigned int const interferingPivotIndex = (*it2);
//                 float const distance13 = _sparseMatrix->getDistance_LocalAndGlobal(queryIndex, interferingPivotIndex, localSparseMatrix);

//                 if (interferingPivotIndex == queryIndex || interferingPivotIndex == activeChildIndex) continue;

//                 // BRUTE FORCE INTERFERENCE CHECK
//                 if (distance13 < distance12 - (2 * queryRadius + activeChildRadius)) {
//                     float const distance23 = _sparseMatrix->getDistance_LocalAndGlobal(activeChildIndex, interferingPivotIndex, localSparseMatrix);
//                     if (distance23 < distance12 - (queryRadius + 2 * activeChildRadius)) {
//                         flag_isNeighbor = false;
//                         break;
//                     }
//                 }
//             }

//             // IF NOT INTERFERED WITH, ADD AS NEIGHBOR
//             if (flag_isNeighbor) {
//                 localFinerPivotLayerNeighborsList.insert(activeChildIndex);
//             }
//         }

// #pragma omp critical(pivotLayerNeighborsListUpdate)
//         { finerLayerPivotNeighborList->insert(localFinerPivotLayerNeighborsList.begin(), localFinerPivotLayerNeighborsList.end()); }

// #pragma omp critical(sparseMatrixUpdate)
//         { _sparseMatrix->updateSparseMatrixWithLocal(localSparseMatrix); }
//     }

//     // sort neighbors to always have them in same order despite threads
//     std::vector<unsigned int> finerPivotLayerNeighborsCopy(finerLayerPivotNeighborList->begin(), finerLayerPivotNeighborList->end());
//     std::sort(finerPivotLayerNeighborsCopy.begin(), finerPivotLayerNeighborsCopy.end());
//     finerLayerPivotNeighborList->clear();
//     finerLayerPivotNeighborList->insert(finerPivotLayerNeighborsCopy.begin(), finerPivotLayerNeighborsCopy.end());

//     return;
// }