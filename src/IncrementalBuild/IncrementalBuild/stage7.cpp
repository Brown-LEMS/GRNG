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

void IB1_collectTopLayerActiveViablePivots(int pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                          RankOrderedList const& rankOrderedTopLayerPivotsList,
                                          tsl::sparse_set<unsigned int>& activeViablePivotsList);
void IB2_collectActiveViablePivots(unsigned int const queryIndex, int const pivotLayerIndex, int const layerIndex,
                                  std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                  tsl::sparse_set<unsigned int> const& activeViableParentsList,
                                  tsl::sparse_set<unsigned int>& activeViableChildrenList);
void IB3_collectActiveViableFinePivots(unsigned int const queryIndex, int const pivotLayerIndex,
                                      std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                      tsl::sparse_set<unsigned int> const& activeViableParentsList,
                                      tsl::sparse_set<unsigned int>& activeViableChildrenList);
void IB4_checkFinePivotsForInterference(unsigned int const queryIndex, int const pivotLayerIndex,
                                       std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                       tsl::sparse_set<unsigned int> const& activeViablePivotsList, LinksMap& invalidatedLinks);
void IB4_checkFinePivotForInterference_multithreading(unsigned int const queryIndex, unsigned int const activePivotIndex, int const pivotLayerIndex,
                                                     std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, LinksMap& localInvalidatedLinks,
                                                     LinksMap& localCheckedLinks, SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix);


/**
 * @brief  Link Validation: Top-down approach to find links invalidated by addition of Query
 *
 * @param queryStruct
 */
void IncrementalBuild::stage7_multithreading(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // initializations
    unsigned int const queryIndex = queryStruct.queryIndex;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    RankOrderedList const& rankOrderedTopLayerPivotsList = queryStruct.rankOrderedTopLayerPivotsList;
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> activeViablePivotsListMap;
    LinksMap invalidatedLinks;

    // top-down approach to collect pivots whose pivot domain may contain exemplars that can be invalidated by query
    // collect list of pivots in Layer pivotLayerIndex-1 who may have children interfered with by Query
    for (int layerIndex = 0; layerIndex < pivotLayerIndex; layerIndex++) {
        activeViablePivotsListMap[layerIndex] = {};

        if (layerIndex == 0) {
            IB1_collectTopLayerActiveViablePivots(pivotLayerIndex, _pivotLayers, rankOrderedTopLayerPivotsList, activeViablePivotsListMap[0]);
        } else {
            IB2_collectActiveViablePivots(queryIndex, pivotLayerIndex, layerIndex, _pivotLayers, activeViablePivotsListMap[layerIndex - 1],
                                         activeViablePivotsListMap[layerIndex]);
        }
    }

    // create list of pivots in Layer pivotLayerIndex that Query may interfere with
    tsl::sparse_set<unsigned int> activeViablePivotsList{};
    activeViablePivotsListMap[pivotLayerIndex] = {};
    IB3_collectActiveViableFinePivots(queryIndex, pivotLayerIndex, _pivotLayers, activeViablePivotsListMap[pivotLayerIndex - 1],
                                     activeViablePivotsList);

    std::vector<unsigned int> activeViablePivotsListVector(activeViablePivotsList.begin(), activeViablePivotsList.end());

    // get links invalidated by the addition of Q in parallel
#pragma omp parallel num_threads(_numThreads)
    {
        LinksMap localInvalidatedLinks;
        LinksMap localCheckedLinks;
        SparseMatrixLocal localSparseMatrix = _sparseMatrix->getLocalCopy();

        // check links of each remaining pivot for invalidation by query
#pragma omp for schedule(static)
        for (int i = 0; i < (int)activeViablePivotsListVector.size(); i++) {
            unsigned int const activePivotIndex = activeViablePivotsListVector[i];
            IB4_checkFinePivotForInterference_multithreading(queryIndex, activePivotIndex, pivotLayerIndex, _pivotLayers, localInvalidatedLinks,
                                                            localCheckedLinks, *_sparseMatrix, localSparseMatrix);
        }

// update with local variants
#pragma omp critical(invalidatedLinksUpdate)
        { invalidatedLinks.merge(localInvalidatedLinks); }

#pragma omp critical(sparseMatrixUpdate)
        { _sparseMatrix->updateSparseMatrixWithLocal(localSparseMatrix); }
    }

    queryStruct.invalidatedPivotLinks = invalidatedLinks;

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[7] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[7] += (dEnd_stage - dStart_stage);
}

/**
 * @brief Check links of activePivotIndex for interference by queryIndex
 *
 * @param queryIndex
 * @param activePivotIndex
 * @param pivotLayerIndex
 * @param pivotLayers
 * @param localInvalidatedLinks
 * @param localCheckedLinks
 * @param sparseMatrix
 * @param localSparseMatrix
 */
void IB4_checkFinePivotForInterference_multithreading(unsigned int const queryIndex, unsigned int const activePivotIndex, int const pivotLayerIndex,
                                                     std::shared_ptr<PivotLayerHierarchy> const& pivotLayers, LinksMap& localInvalidatedLinks,
                                                     LinksMap& localCheckedLinks, SparseMatrix& sparseMatrix, SparseMatrixLocal& localSparseMatrix) {
    float const activePivotRadius = (*pivotLayers)[pivotLayerIndex].get_pivotRadius(activePivotIndex);

    // check query against pivot's max link distance if possible
    std::pair<bool, float> const distance13_pair =
        sparseMatrix.getDistanceIfAvailable_LocalAndGlobal(activePivotIndex, queryIndex, localSparseMatrix);
    if (distance13_pair.first) {
        IndexPairDistanceStruct const& umax_activePivotIndex = (*pivotLayers)[pivotLayerIndex].get_umaxMaxLinkDistance(activePivotIndex);
        if (!umax_activePivotIndex.isInvalid()) {
            if (distance13_pair.second >= umax_activePivotIndex.distance) {
                return;
            }
        }
    }

    // check individual links by neighbors of activePivotIndex
    tsl::sparse_set<unsigned int> const& activePivotNeighborsList = (*pivotLayers)[pivotLayerIndex].get_pivotLayerNeighbors(activePivotIndex);
    tsl::sparse_set<unsigned int>::const_iterator it2;
    for (it2 = activePivotNeighborsList.begin(); it2 != activePivotNeighborsList.end(); it2++) {
        unsigned int const activePivotNeighborIndex = (*it2);

        // ensure link hasn't already been checked by the other side
        if (localCheckedLinks.checkAndAddLink(activePivotIndex, activePivotNeighborIndex)) continue;

        float const activePivotNeighborRadius = (*pivotLayers)[pivotLayerIndex].get_pivotRadius(activePivotNeighborIndex);
        float const distance12 = sparseMatrix.getDistance_LocalAndGlobal(activePivotIndex, activePivotNeighborIndex, localSparseMatrix);
        float const activePivotIndexLinkDistance = pivotLayers->_luneRadius(distance12, activePivotRadius, activePivotNeighborRadius);
        float const activePivotNeighborLinkDistance = pivotLayers->_luneRadius(distance12, activePivotNeighborRadius, activePivotRadius);

        // check if link has a lune -> if not, cannot interfere
        if (activePivotIndexLinkDistance <= 0 || activePivotNeighborLinkDistance <= 0) continue;

        // lastly, brute force check for interference
        float const distance13 = sparseMatrix.getDistance_LocalAndGlobal(activePivotIndex, queryIndex, localSparseMatrix);
        if (distance13 < activePivotIndexLinkDistance) {
            float const distance23 = sparseMatrix.getDistance_LocalAndGlobal(activePivotNeighborIndex, queryIndex, localSparseMatrix);
            if (distance23 < activePivotNeighborLinkDistance) {
                localInvalidatedLinks.checkAndAddLink(activePivotIndex, activePivotNeighborIndex);  // link invalidated by Q!
            }
        }
    }

    return;
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

/**
 * @brief  Link Validation: Top-down approach to find links invalidated by addition of Query
 *
 * @param queryStruct
 */
void IncrementalBuild::stage7(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // initializations
    unsigned int const queryIndex = queryStruct.queryIndex;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    RankOrderedList const& rankOrderedTopLayerPivotsList = queryStruct.rankOrderedTopLayerPivotsList;
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> activeViablePivotsListMap;
    LinksMap invalidatedLinks;

    // top-down approach to collect pivots whose pivot domain may contain exemplars that can be invalidated by query
    // collect list of pivots in Layer pivotLayerIndex-1 who may have children interfered with by Query
    for (int layerIndex = 0; layerIndex < pivotLayerIndex; layerIndex++) {
        activeViablePivotsListMap[layerIndex] = {};

        if (layerIndex == 0) {
            IB1_collectTopLayerActiveViablePivots(pivotLayerIndex, _pivotLayers, rankOrderedTopLayerPivotsList, activeViablePivotsListMap[0]);
        } else {
            IB2_collectActiveViablePivots(queryIndex, pivotLayerIndex, layerIndex, _pivotLayers, activeViablePivotsListMap[layerIndex - 1],
                                         activeViablePivotsListMap[layerIndex]);
        }
    }

    // create list of pivots in Layer pivotLayerIndex that Query may interfere with
    activeViablePivotsListMap[pivotLayerIndex] = {};
    IB3_collectActiveViableFinePivots(queryIndex, pivotLayerIndex, _pivotLayers, activeViablePivotsListMap[pivotLayerIndex - 1],
                                     activeViablePivotsListMap[pivotLayerIndex]);

    // get links invalidated by the addition of Q
    IB4_checkFinePivotsForInterference(queryIndex, pivotLayerIndex, _pivotLayers, activeViablePivotsListMap[pivotLayerIndex], invalidatedLinks);
    queryStruct.invalidatedPivotLinks = invalidatedLinks;

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[7] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[7] += (dEnd_stage - dStart_stage);
}



/**
 * @brief Collect all top layer pivots that have children that query can interfere with.
 *
 * @param pivotLayerIndex layer that query will be added to
 * @param pivotLayers hierarchy data structure
 * @param rankOrderedTopLayerPivotsList ranked list of top layer pivots from Stage 0
 * @param activeViablePivotsList all active top layer pivots
 */
void IB1_collectTopLayerActiveViablePivots(int pivotLayerIndex, std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                          RankOrderedList const& rankOrderedTopLayerPivotsList,
                                          tsl::sparse_set<unsigned int>& activeViablePivotsList) {
    activeViablePivotsList.clear();

    // umaxMaximumOfDescendantMaxLinkDistance to find furthest reach of query
    IndexPairDistanceStruct umax_layer = (*pivotLayers)[0].get_umaxMaximumOfDescendantMaxLinkDistance(pivotLayerIndex);
    if (umax_layer.isInvalid()) return;  // no links to break in layer

    // consider top layer pivots by ranked list
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = rankOrderedTopLayerPivotsList.list.begin(); it1 != rankOrderedTopLayerPivotsList.list.end(); it1++) {
        unsigned int const topLayerPivot = (*it1).index;
        float const distance = (*it1).distance;
        bool flag_add = true;

        // query out of range of all links in layer, can break since ranked list
        if (distance >= umax_layer.distance) break;

        // check if query out of range of this pivot
        IndexPairDistanceStruct const& umax_topLayerPivot = (*pivotLayers)[0].get_umaxDescendantMaxLinkDistance(topLayerPivot, pivotLayerIndex);
        if (!umax_topLayerPivot.isInvalid()) {
            if (distance >= umax_topLayerPivot.distance) {
                flag_add = false;
            }
        }

        // add if pivot passes umax tests, has children within reach of query
        if (flag_add) activeViablePivotsList.insert(topLayerPivot);
    }

    return;
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
void IB2_collectActiveViablePivots(unsigned int const queryIndex, int const pivotLayerIndex, int const layerIndex,
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

    // check all children for passing umax tests
    tsl::sparse_set<unsigned int>::const_iterator it4;
    for (it4 = candidateActiveViablePivotsList.begin(); it4 != candidateActiveViablePivotsList.end(); it4++) {
        unsigned int candidateActiveChildIndex = (*it4);
        float const distance = pivotLayers->getDistance(queryIndex, candidateActiveChildIndex);
        bool flag_add = true;

        // check against umax
        IndexPairDistanceStruct const& umax_activeChildIndex =
            (*pivotLayers)[layerIndex].get_umaxDescendantMaxLinkDistance(candidateActiveChildIndex, pivotLayerIndex);
        if (!umax_activeChildIndex.isInvalid()) {
            if (distance >= umax_activeChildIndex.distance) {
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
void IB3_collectActiveViableFinePivots(unsigned int const queryIndex, int const pivotLayerIndex,
                                      std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                      tsl::sparse_set<unsigned int> const& activeViableParentsList,
                                      tsl::sparse_set<unsigned int>& activeViableChildrenList) {
    activeViableChildrenList.clear();

    // create list of children of active pivots in the coarse layer
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeViableParentsList.begin(); it1 != activeViableParentsList.end(); it1++) {
        unsigned int const activeParentIndex = (*it1);

        // consider all children of parent
        tsl::sparse_set<unsigned int> const& activeParentIndexChildren =
            (*pivotLayers)[pivotLayerIndex - 1].get_descendantPivotIndices(activeParentIndex, pivotLayerIndex);

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activeParentIndexChildren.begin(); it2 != activeParentIndexChildren.end(); it2++) {
            unsigned int const childIndex = (*it2);
            bool flag_add = true;

            // check that all parents of childIndex are in the activeViable pivots list
            tsl::sparse_set<unsigned int> const& childIndexParents =
                (*pivotLayers)[pivotLayerIndex].get_ancestorPivotIndices(childIndex, pivotLayerIndex - 1);

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
void IB4_checkFinePivotsForInterference(unsigned int const queryIndex, int const pivotLayerIndex,
                                       std::shared_ptr<PivotLayerHierarchy> const& pivotLayers,
                                       tsl::sparse_set<unsigned int> const& activeViablePivotsList, LinksMap& invalidatedLinks) {
    invalidatedLinks.clear();
    LinksMap checkedLinks;

    // check each active pivot for link interference by Query
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeViablePivotsList.begin(); it1 != activeViablePivotsList.end(); it1++) {
        unsigned int const activePivotIndex = (*it1);
        float const activePivotRadius = (*pivotLayers)[pivotLayerIndex].get_pivotRadius(activePivotIndex);

        // check query against pivot's max link distance if possible
        std::pair<bool, float> const distance13_pair = pivotLayers->getDistanceIfAvailable(activePivotIndex, queryIndex);
        if (distance13_pair.first) {
            IndexPairDistanceStruct const& umax_activePivotIndex = (*pivotLayers)[pivotLayerIndex].get_umaxMaxLinkDistance(activePivotIndex);
            if (!umax_activePivotIndex.isInvalid()) {
                if (distance13_pair.second >= umax_activePivotIndex.distance) {
                    continue;
                }
            }
        }

        // check individual links by neighbors of activePivotIndex
        tsl::sparse_set<unsigned int> const& activePivotNeighborsList = (*pivotLayers)[pivotLayerIndex].get_pivotLayerNeighbors(activePivotIndex);
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = activePivotNeighborsList.begin(); it2 != activePivotNeighborsList.end(); it2++) {
            unsigned int const activePivotNeighborIndex = (*it2);

            // ensure link hasn't already been checked by the other side
            if (checkedLinks.checkAndAddLink(activePivotIndex, activePivotNeighborIndex)) continue;

            float const activePivotNeighborRadius = (*pivotLayers)[pivotLayerIndex].get_pivotRadius(activePivotNeighborIndex);
            float const distance12 = pivotLayers->getDistance(activePivotIndex, activePivotNeighborIndex);
            float const activePivotIndexLinkDistance = pivotLayers->_luneRadius(distance12, activePivotRadius, activePivotNeighborRadius);
            float const activePivotNeighborLinkDistance = pivotLayers->_luneRadius(distance12, activePivotNeighborRadius, activePivotRadius);

            // check if link has a lune -> if not, cannot interfere
            if (activePivotIndexLinkDistance <= 0 || activePivotNeighborLinkDistance <= 0) continue;

            // lastly, brute force check for interference
            float const distance13 = pivotLayers->getDistance(activePivotIndex, queryIndex);
            if (distance13 < activePivotIndexLinkDistance) {
                float const distance23 = pivotLayers->getDistance(activePivotNeighborIndex, queryIndex);
                if (distance23 < activePivotNeighborLinkDistance) {
                    invalidatedLinks.checkAndAddLink(activePivotIndex, activePivotNeighborIndex);  // link invalidated by Q!
                }
            }
        }
    }

    return;
}



/**
 * @brief Brute force approach to link invalidation. Check each link to be invalidated by Q. Good for debugging.
 *
 * @param QueryStruct
 */
void IncrementalBuild::stage7_bruteForce(QueryStructIncrementalBuild& QueryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // GET INFORMATION FROM QUERY OBJECT
    unsigned int const queryIndex = QueryStruct.queryIndex;
    int const pivotLayerIndex = QueryStruct.pivotLayerIndex;
    LinksMap invalidatedLinks;

    // CHECK LINKS OF EACH PIVOT FOR INTERFERENCE
    tsl::sparse_set<unsigned int> const& activePivotsList = *((*_pivotLayers)[pivotLayerIndex].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activePivotsList.begin(); it1 != activePivotsList.end(); it1++) {
        unsigned int const activePivotIndex = (*it1);
        float const activePivotRadius = (*_pivotLayers)[pivotLayerIndex].get_pivotRadius(activePivotIndex);
        float const distance13 = _pivotLayers->getDistance(queryIndex, activePivotIndex);

        // CHECK ALL LINKS IF INTERFERED WITH
        tsl::sparse_set<unsigned int> const& neighborsList = (*_pivotLayers)[pivotLayerIndex].get_pivotLayerNeighbors(activePivotIndex);
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = neighborsList.begin(); it2 != neighborsList.end(); it2++) {
            // neighbor of active pivot information
            unsigned int const activePivotNeighborIndex = (*it2);
            float const activePivotNeighborRadius = (*_pivotLayers)[pivotLayerIndex].get_pivotRadius(activePivotNeighborIndex);
            float const distance12 = _pivotLayers->getDistance(activePivotIndex, activePivotNeighborIndex);

            // IF NO LUNE, CANNOT BE INTERFERED WITH
            if (distance12 - (2 * activePivotRadius + activePivotNeighborRadius) <= 0 ||
                distance12 - (activePivotRadius + 2 * activePivotNeighborRadius) <= 0) {
                continue;
            }

            // BRUTE FORCE CHECK FOR QUERY INTERFERENCE
            if (distance13 < distance12 - (2 * activePivotRadius + activePivotNeighborRadius)) {
                float const distance23 = _pivotLayers->getDistance(activePivotNeighborIndex, queryIndex);

                if (distance23 < distance12 - (activePivotRadius + 2 * activePivotNeighborRadius)) {
                    invalidatedLinks.checkAndAddLink(activePivotIndex, activePivotNeighborIndex);
                }
            }
        }
    }

    // UPDATE QUERY OBJECT INFORMATION
    QueryStruct.invalidatedPivotLinks = invalidatedLinks;

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[7] = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[7] = (dEnd_stage - dStart_stage);

    return;
}