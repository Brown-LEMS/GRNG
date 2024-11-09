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

namespace IncrementalBuild_Stage9 {

void collectTopLayerActiveViablePivots(int pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                       RankOrderedList const& rankOrderedTopLayerPivotsList, tsl::sparse_set<unsigned int>& activeViablePivotsList);

void collectActiveViablePivots(unsigned int const queryIndex, int const pivotLayerIndex, int const layerIndex,
                               std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, tsl::sparse_set<unsigned int> const& activeViableParentsList,
                               tsl::sparse_set<unsigned int>& activeViableChildrenList);

void collectActiveViableFinePivots(unsigned int const queryIndex, int const pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                   tsl::sparse_set<unsigned int> const& activeViableParentsList,
                                   tsl::sparse_set<unsigned int>& activeViableChildrenList);

void checkFinePivotsForInterference(unsigned int const queryIndex, int const pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                    tsl::sparse_set<unsigned int> const& activeViablePivotsList, LinksMap& invalidatedLinks);

void checkFinePivotForInterference_multithreading(unsigned int const queryIndex, unsigned int const activePivotIndex, int const pivotLayerIndex,
                                                  std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, LinksMap& localInvalidatedLinks,
                                                  LinksMap& localCheckedLinks, SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix);

bool vmaxLinkCheck(unsigned int const queryIndex, unsigned int const finePivotIndex, float const distance, int const finerLayerIndex,
                   std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);

};  // namespace IncrementalBuild_Stage9

void IncrementalBuild::stage9(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // GET ALL STRUCTURES FROM QIB
    unsigned int const queryIndex = queryStruct.queryIndex;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    int const finerLayerIndex = pivotLayerIndex + 1;
    RankOrderedList const& rankOrderedTopLayerPivotsList = queryStruct.rankOrderedTopLayerPivotsList;
    LinksMap& invalidatedFineCoarseNeighborsLinks = queryStruct.invalidatedFineCoarseNeighborsLinks;

    // need to remove Q from list of rank ordered pivots
    // RankOrderedList rankOrderedTopLayerPivots;
    // if (pivotLayerIndex == 0) {
    //     std::vector<IndexDistanceStruct>::const_iterator it03;
    //     for (it03 = rankOrderedTopLayerPivotsList.list.begin(); it03 != rankOrderedTopLayerPivotsList.list.end(); it03++) {
    //         unsigned int const topIndex = (*it03).index;
    //         if (topIndex != queryIndex) {
    //             rankOrderedTopLayerPivots.add((*it03).index,(*it03).distance);
    //         }
    //     }
    //     rankOrderedTopLayerPivots.sort();
    // } else {
    //     rankOrderedTopLayerPivots = rankOrderedTopLayerPivotsList;
    // }

    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> activeViablePivotsListMap{};
    tsl::sparse_set<unsigned int> finePivotsToCheck{};
    for (int layerIndex = 0; layerIndex < finerLayerIndex; layerIndex++) {
        activeViablePivotsListMap[layerIndex] = {};

        if (layerIndex == 0) {
            IncrementalBuild_Stage9::collectTopLayerActiveViablePivots(finerLayerIndex, _pivotLayers, rankOrderedTopLayerPivotsList,
                                                                       activeViablePivotsListMap[0]);
        } else {
            IncrementalBuild_Stage9::collectActiveViablePivots(queryIndex, finerLayerIndex, layerIndex, _pivotLayers,
                                                               activeViablePivotsListMap[layerIndex - 1], activeViablePivotsListMap[layerIndex]);
        }
    }

    // create list of pivots in Layer pivotLayerIndex that Query may interfere with
    tsl::sparse_set<unsigned int> activeViablePivotsList{};
    IncrementalBuild_Stage9::collectActiveViableFinePivots(queryIndex, finerLayerIndex, _pivotLayers, activeViablePivotsListMap[pivotLayerIndex],
                                                           activeViablePivotsList);

    IncrementalBuild_Stage9::checkFinePivotsForInterference(queryIndex, finerLayerIndex, _pivotLayers, activeViablePivotsList,
                                                            invalidatedFineCoarseNeighborsLinks);

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[9] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[9] += (dEnd_stage - dStart_stage);
}

void IncrementalBuild_Stage9::collectTopLayerActiveViablePivots(int finerLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                                RankOrderedList const& rankOrderedTopLayerPivotsList,
                                                                tsl::sparse_set<unsigned int>& activeViablePivotsList) {
    activeViablePivotsList.clear();

    // umaxMaximumOfDescendantMaxLinkDistance to find furthest reach of query
    IndexPairDistanceStruct vmax_layer = (*pivotLayers)[0].get_vmaxMaximumOfDescendantMaxCoarseLinkDistance(finerLayerIndex);
    if (vmax_layer.isInvalid()) return;  // no links to break in layer

    // consider top layer pivots by ranked list
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = rankOrderedTopLayerPivotsList.list.begin(); it1 != rankOrderedTopLayerPivotsList.list.end(); it1++) {
        unsigned int const topLayerPivot = (*it1).index;
        float const distance = (*it1).distance;
        bool flag_add = true;

        // query out of range of all links in layer, can break since ranked list
        if (distance >= vmax_layer.distance) break;

        // check if query out of range of this pivot
        IndexPairDistanceStruct const& vmax_topLayerPivot = (*pivotLayers)[0].get_vmaxDescendantMaxCoarseLinkDistance(topLayerPivot, finerLayerIndex);
        if (!vmax_topLayerPivot.isInvalid()) {
            if (distance >= vmax_topLayerPivot.distance) {
                flag_add = false;
            }
        }

        // add if pivot passes umax tests, has children within reach of query
        if (flag_add) activeViablePivotsList.insert(topLayerPivot);
    }
}

/**
 * @brief Collect pivots of lower layers that have children that query may interfere with
 *
 * @param queryIndex
 * @param pivotLayerIndex layer query is being added to
 * @param layerIndex layer of pivot children we are collecting
 * @param pivotLayers
 * @param activeViableParentsList list of active parents
 * @param activeViableChildrenList list fo pivot children we are collecting
 */
void IncrementalBuild_Stage9::collectActiveViablePivots(unsigned int const queryIndex, int const finerLayerIndex, int const layerIndex,
                                                        std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                        tsl::sparse_set<unsigned int> const& activeViableParentsList,
                                                        tsl::sparse_set<unsigned int>& activeViableChildrenList) {
    activeViableChildrenList.clear();

    // create list of children of active pivots in the above layer
    tsl::sparse_set<unsigned int> candidateActiveViablePivotsList{};
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeViableParentsList.begin(); it1 != activeViableParentsList.end(); it1++) {
        unsigned int const activeParentIndex = (*it1);

        // consider all children of parent
        tsl::sparse_set<unsigned int> const& activeParentIndexChildren =
            (*pivotLayers)[layerIndex - 1].get_descendantPivotIndices(activeParentIndex, layerIndex);
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activeParentIndexChildren.begin(); it2 != activeParentIndexChildren.end(); it2++) {
            unsigned int const childIndex = (*it2);
            bool flag_add = true;

            // check that all parents of childIndex are in the activeViable pivots list
            tsl::sparse_set<unsigned int> const& childIndexParents = (*pivotLayers)[layerIndex].get_ancestorPivotIndices(childIndex, layerIndex - 1);
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = childIndexParents.begin(); it3 != childIndexParents.end(); it3++) {
                unsigned int const childParentIndex = (*it3);

                if (activeViableParentsList.find(childParentIndex) == activeViableParentsList.end()) {
                    flag_add = false;
                    break;
                }
            }

            // add child as potentially active pivot
            if (flag_add) candidateActiveViablePivotsList.insert(childIndex);
        }
    }

    // check all children for passing vmax tests
    tsl::sparse_set<unsigned int>::const_iterator it4;
    for (it4 = candidateActiveViablePivotsList.begin(); it4 != candidateActiveViablePivotsList.end(); it4++) {
        unsigned int candidateActiveChildIndex = (*it4);
        float const distance = pivotLayers->getDistance(queryIndex, candidateActiveChildIndex);
        bool flag_add = true;

        // check against umax
        IndexPairDistanceStruct const& vmax_activeChildIndex =
            (*pivotLayers)[layerIndex].get_vmaxDescendantMaxCoarseLinkDistance(candidateActiveChildIndex, finerLayerIndex);
        if (!vmax_activeChildIndex.isInvalid()) {
            if (distance >= vmax_activeChildIndex.distance) {
                flag_add = false;
            }
        }

        if (flag_add) activeViableChildrenList.insert(candidateActiveChildIndex);
    }

    return;
}

/**
 * @brief given list of coarse pivots, collect pivots in query's layer that it may interfere with
 *
 * @param queryIndex
 * @param pivotLayerIndex // layer query is being added to
 * @param pivotLayers
 * @param activeViableParentsList // list of coarse active pivots whose domain contains pivots Q may interfere with
 * @param activeViableChildrenList // list of pivots in query's layer it may interfere with
 */
void IncrementalBuild_Stage9::collectActiveViableFinePivots(unsigned int const queryIndex, int const finerLayerIndex,
                                                            std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                            tsl::sparse_set<unsigned int> const& activeViableParentsList,
                                                            tsl::sparse_set<unsigned int>& activeViableChildrenList) {
    activeViableChildrenList.clear();

    // create list of children of active pivots in the coarse layer
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeViableParentsList.begin(); it1 != activeViableParentsList.end(); it1++) {
        unsigned int const activeParentIndex = (*it1);

        tsl::sparse_set<unsigned int> const& indices = *(*pivotLayers)[finerLayerIndex - 1].get_pivotIndices_ptr();
        // consider all children of parent
        tsl::sparse_set<unsigned int> const& activeParentIndexChildren =
            (*pivotLayers)[finerLayerIndex - 1].get_descendantPivotIndices(activeParentIndex, finerLayerIndex);

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activeParentIndexChildren.begin(); it2 != activeParentIndexChildren.end(); it2++) {
            unsigned int const childIndex = (*it2);
            bool flag_add = true;

            // check that all parents of childIndex are in the activeViable pivots list
            tsl::sparse_set<unsigned int> const& childIndexParents =
                (*pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(childIndex, finerLayerIndex - 1);

            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = childIndexParents.begin(); it3 != childIndexParents.end(); it3++) {
                unsigned int const childParentIndex = (*it3);

                if (activeViableParentsList.find(childParentIndex) == activeViableParentsList.end()) {
                    flag_add = false;
                    break;
                }
            }

            // add child as potentially active pivot
            if (flag_add) activeViableChildrenList.insert(childIndex);
        }
    }

    return;
}

/**
 * @brief gather all invalidated links by query based on pivots it may interfere with
 *
 * @param queryIndex
 * @param pivotLayerIndex
 * @param pivotLayers
 * @param activeViablePivotsList list of pivots whose links Q might interfere with
 * @param invalidatedLinks list of invalidated links
 */
void IncrementalBuild_Stage9::checkFinePivotsForInterference(unsigned int const queryIndex, int const finerLayerIndex,
                                                             std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                             tsl::sparse_set<unsigned int> const& activeViablePivotsList,
                                                             LinksMap& invalidatedLinks) {
    invalidatedLinks.clear();

    // check each active FINE pivot for coarse link interference by Query
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeViablePivotsList.begin(); it1 != activeViablePivotsList.end(); it1++) {
        unsigned int const activePivotIndex = (*it1);
        float const activePivotRadius = (*pivotLayers)[finerLayerIndex].get_pivotRadius(activePivotIndex);
        float const distance13 = (*pivotLayers).getDistance(activePivotIndex, queryIndex);

        // vMax check to see if we need to check all coarse neighbors of finePivot
        if (!IncrementalBuild_Stage9::vmaxLinkCheck(queryIndex, activePivotIndex, distance13, finerLayerIndex, pivotLayers)) {
            continue;
        }

        // check individual links by neighbors of activePivotIndex
        tsl::sparse_set<unsigned int> const& activeCoarseNeighborsList =
            (*pivotLayers)[finerLayerIndex].get_coarsePivotLayerNeighbors(activePivotIndex);
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activeCoarseNeighborsList.begin(); it2 != activeCoarseNeighborsList.end(); it2++) {
            unsigned int const activePivotNeighborIndex = (*it2);

            float const activePivotNeighborRadius = (*pivotLayers)[finerLayerIndex - 1].get_pivotRadius(activePivotNeighborIndex);
            float const distance12 = pivotLayers->getDistance(activePivotIndex, activePivotNeighborIndex);
            float const activePivotIndexLinkDistance = pivotLayers->_luneRadius(distance12, activePivotRadius, activePivotNeighborRadius);
            float const activePivotNeighborLinkDistance = pivotLayers->_luneRadius(distance12, activePivotNeighborRadius, activePivotRadius);

            // check if link has a lune -> if not, cannot interfere
            if (activePivotIndexLinkDistance <= 0 || activePivotNeighborLinkDistance <= 0) continue;

            // lastly, brute force check for interference
            if (distance13 < activePivotIndexLinkDistance) {
                float const distance23 = pivotLayers->getDistance(activePivotNeighborIndex, queryIndex);
                if (distance23 < activePivotNeighborLinkDistance) {
                    invalidatedLinks.addOrderedIndexLink(activePivotIndex, activePivotNeighborIndex);  // coarse link invalidated by Q!
                }
            }
        }
    }

    return;
}

bool IncrementalBuild_Stage9::vmaxLinkCheck(unsigned int const queryIndex, unsigned int const finePivotIndex, float const distance,
                                            int const finerLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers) {
    bool flag_inReach = true;

    IndexPairDistanceStruct const& vmax = (*pivotLayers)[finerLayerIndex].get_vmaxMaxCoarseLinkDistance(finePivotIndex);

    if (!vmax.isInvalid()) {
        if (distance > vmax.distance) {
            flag_inReach = false;
        }
    }

    return flag_inReach;
}

/**
 *
 *
 *
 *
 *
 *
 *
 *                      MULTITHREADING
 *
 *
 *
 *
 *
 *
 *
 */

namespace IncrementalBuild_Stage9_Multithreading {

void checkFinePivotsForInterference(unsigned int const queryIndex, unsigned int const activePivotIndex, int const finerLayerIndex,
                                    std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, LinksMap& invalidatedLinks, SparseMatrix& sparseMatrix,
                                    SparseMatrixLocal& localSparseMatrix);

};  // namespace IncrementalBuild_Stage9_Multithreading

/**
 * @brief Stage 9: Find Coarse-Fine Links Invalidated by addition of Q
 *
 * @param queryStruct
 */
void IncrementalBuild::stage9_multithreading(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // GET ALL STRUCTURES FROM QIB
    unsigned int const queryIndex = queryStruct.queryIndex;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    int const finerLayerIndex = pivotLayerIndex + 1;
    RankOrderedList const& rankOrderedTopLayerPivotsList = queryStruct.rankOrderedTopLayerPivotsList;
    LinksMap& invalidatedFineCoarseNeighborsLinks = queryStruct.invalidatedFineCoarseNeighborsLinks;

    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> activeViablePivotsListMap{};
    tsl::sparse_set<unsigned int> finePivotsToCheck{};
    for (int layerIndex = 0; layerIndex < finerLayerIndex; layerIndex++) {
        activeViablePivotsListMap[layerIndex] = {};

        if (layerIndex == 0) {
            IncrementalBuild_Stage9::collectTopLayerActiveViablePivots(finerLayerIndex, _pivotLayers, rankOrderedTopLayerPivotsList,
                                                                       activeViablePivotsListMap[0]);
        } else {
            IncrementalBuild_Stage9::collectActiveViablePivots(queryIndex, finerLayerIndex, layerIndex, _pivotLayers,
                                                               activeViablePivotsListMap[layerIndex - 1], activeViablePivotsListMap[layerIndex]);
        }
    }

    // create list of pivots in Layer pivotLayerIndex that Query may interfere with
    tsl::sparse_set<unsigned int> activeViablePivotsList{};
    IncrementalBuild_Stage9::collectActiveViableFinePivots(queryIndex, finerLayerIndex, _pivotLayers, activeViablePivotsListMap[pivotLayerIndex],
                                                           activeViablePivotsList);

    std::vector<unsigned int> activeViablePivotsListVector(activeViablePivotsList.begin(), activeViablePivotsList.end());
#pragma omp parallel num_threads(_numThreads)
    {
        LinksMap localInvalidatedLinks;
        SparseMatrixLocal localSparseMatrix = _sparseMatrix->getLocalCopy();

        // check links of each remaining pivot for invalidation by query
#pragma omp for schedule(static)
        for (int i = 0; i < (int)activeViablePivotsListVector.size(); i++) {
            unsigned int const activePivotIndex = activeViablePivotsListVector[i];
            IncrementalBuild_Stage9_Multithreading::checkFinePivotsForInterference(queryIndex, activePivotIndex, finerLayerIndex, _pivotLayers,
                                                                                   localInvalidatedLinks, *_sparseMatrix, localSparseMatrix);
        }

// update with local variants
#pragma omp critical(invalidatedLinksUpdate)
        { invalidatedFineCoarseNeighborsLinks.merge(localInvalidatedLinks); }

#pragma omp critical(sparseMatrixUpdate)
        { _sparseMatrix->updateSparseMatrixWithLocal(localSparseMatrix); }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[9] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[9] += (dEnd_stage - dStart_stage);
    return;
}

/**
 * @brief Check all coarse links of activePivotIndex for interference by Q. uses multithreading
 *
 * @param queryIndex
 * @param activePivotIndex
 * @param finerLayerIndex
 * @param pivotLayers
 * @param invalidatedLinks
 * @param sparseMatrix
 * @param localSparseMatrix
 */
void IncrementalBuild_Stage9_Multithreading::checkFinePivotsForInterference(unsigned int const queryIndex, unsigned int const activePivotIndex,
                                                                            int const finerLayerIndex,
                                                                            std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                                                            LinksMap& invalidatedLinks, SparseMatrix& sparseMatrix,
                                                                            SparseMatrixLocal& localSparseMatrix) {
    float const activePivotRadius = (*pivotLayers)[finerLayerIndex].get_pivotRadius(activePivotIndex);
    float const distance13 = sparseMatrix.getDistance_LocalAndGlobal(activePivotIndex, queryIndex, localSparseMatrix);

    // vMax check to see if we need to check all coarse neighbors of finePivot
    if (!IncrementalBuild_Stage9::vmaxLinkCheck(queryIndex, activePivotIndex, distance13, finerLayerIndex, pivotLayers)) {
        return;
    }

    // check individual links by neighbors of activePivotIndex
    tsl::sparse_set<unsigned int> const& activeCoarseNeighborsList = (*pivotLayers)[finerLayerIndex].get_coarsePivotLayerNeighbors(activePivotIndex);
    tsl::sparse_set<unsigned int>::const_iterator it2;
    for (it2 = activeCoarseNeighborsList.begin(); it2 != activeCoarseNeighborsList.end(); it2++) {
        unsigned int const activePivotNeighborIndex = (*it2);

        float const activePivotNeighborRadius = (*pivotLayers)[finerLayerIndex - 1].get_pivotRadius(activePivotNeighborIndex);
        float const distance12 = sparseMatrix.getDistance_LocalAndGlobal(activePivotIndex, activePivotNeighborIndex, localSparseMatrix);
        float const activePivotIndexLinkDistance = pivotLayers->_luneRadius(distance12, activePivotRadius, activePivotNeighborRadius);
        float const activePivotNeighborLinkDistance = pivotLayers->_luneRadius(distance12, activePivotNeighborRadius, activePivotRadius);

        // check if link has a lune -> if not, cannot interfere
        if (activePivotIndexLinkDistance <= 0 || activePivotNeighborLinkDistance <= 0) continue;

        // lastly, brute force check for interference
        if (distance13 < activePivotIndexLinkDistance) {
            float const distance23 = sparseMatrix.getDistance_LocalAndGlobal(activePivotNeighborIndex, queryIndex, localSparseMatrix);
            if (distance23 < activePivotNeighborLinkDistance) {
                invalidatedLinks.addOrderedIndexLink(activePivotIndex, activePivotNeighborIndex);  // coarse link invalidated by Q!
            }
        }
    }

    return;
}

/**
 *
 *
 *
 *
 *
 *
 *
 *                      BRUTE FORCE
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
void IncrementalBuild::stage9_bruteForce(QueryStructIncrementalBuild& queryStruct) {
    // GET ALL STRUCTURES FROM QIB
    unsigned int const queryIndex = queryStruct.queryIndex;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    int const finerLayerIndex = pivotLayerIndex + 1;
    LinksMap& invalidatedFineCoarseNeighborsLinks = queryStruct.invalidatedFineCoarseNeighborsLinks;

    // need to check all fine pivots for interference
    tsl::sparse_set<unsigned int> const& finePivotIndices = *(*_pivotLayers)[finerLayerIndex].get_pivotIndices_ptr();
    tsl::sparse_set<unsigned int>::const_iterator it1, it2;
    for (it1 = finePivotIndices.begin(); it1 != finePivotIndices.end(); it1++) {
        unsigned int const finePivotIndex = (*it1);
        float const finePivotRadius = (*_pivotLayers)[finerLayerIndex].get_pivotRadius(finePivotIndex);

        tsl::sparse_set<unsigned int> const& finePivotCoarseNeighbors =
            (*_pivotLayers)[finerLayerIndex].get_coarsePivotLayerNeighbors(finePivotIndex);

        for (it2 = finePivotCoarseNeighbors.begin(); it2 != finePivotCoarseNeighbors.end(); it2++) {
            unsigned int const coarsePivotIndex = (*it2);
            float const coarsePivotRadius = (*_pivotLayers)[pivotLayerIndex].get_pivotRadius(coarsePivotIndex);
            float const distance12 = _sparseMatrix->getDistance(finePivotIndex, coarsePivotIndex);
            bool flag_neighbor = true;

            // check if query interferees!
            float const distance13 = _sparseMatrix->getDistance(finePivotIndex, queryIndex);
            if (distance13 < _luneRadius(distance12, finePivotRadius, coarsePivotRadius)) {  // reverse for check. coarse radius bigger.
                float const distance23 = _sparseMatrix->getDistance(coarsePivotIndex, queryIndex);
                if (distance23 < _luneRadius(distance12, coarsePivotRadius, finePivotRadius)) {
                    flag_neighbor = false;
                }
            }

            // coarse neighbor of finePivotIndex is invalidated!
            if (!flag_neighbor) {
                invalidatedFineCoarseNeighborsLinks.addOrderedIndexLink(finePivotIndex, coarsePivotIndex);  // store ordered
            };
        }
    }

    return;
}
