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
// 2022-06-01
#include <omp.h>

#include "online-search.hpp"

/**
 * @brief Pivot-Exemplar Interactions: check that exemplars are coarse neighbors with all parents of Q.
 *
 * @param queryStruct
 */
void OnlineSearch::stage3_multithreading(QueryStructOnlineSearch& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // initializations
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int> const& coarserPivotLayerNeighborsList = queryStruct.pivotLayerNeighborsListMap[coarserLayerIndex];
    tsl::sparse_set<unsigned int> const& parentsList = queryStruct.ancestorIndicesMap[coarserLayerIndex];
    tsl::sparse_set<unsigned int>& activeChildrenList = queryStruct.activeChildrenList;
    tsl::sparse_set<unsigned int>& finerPivotLayerNeighborsList = queryStruct.pivotLayerNeighborsListMap[finerLayerIndex];

    std::vector<unsigned int> coarserPivotLayerNeighborsVector(coarserPivotLayerNeighborsList.begin(), coarserPivotLayerNeighborsList.end());
    tsl::sparse_set<unsigned int> potentialFineLayerPivotNeighborsList{};
#pragma omp parallel num_threads(_numThreads)
    {
        tsl::sparse_set<unsigned int> localPotentialNeighborsList;

        // loop through each coarse neighbor
#pragma omp for schedule(static)
        for (int it1 = 0; it1 < (int)coarserPivotLayerNeighborsVector.size(); it1++) {
            unsigned int const coarserNeighborIndex = coarserPivotLayerNeighborsVector[it1];

            // consider each child of the coarse neighbor
            tsl::sparse_set<unsigned int> const& coarseNeighborChildrenList =
                (*_pivotLayers)[coarserLayerIndex].get_descendantPivotIndices(coarserNeighborIndex, finerLayerIndex);
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = coarseNeighborChildrenList.begin(); it2 != coarseNeighborChildrenList.end(); it2++) {
                unsigned int const childIndex = (*it2);
                bool flag_addChild = true;

                // check that query is coarse neighbors to all parents of childIndex
                tsl::sparse_set<unsigned int> const& childIndexParentList =
                    (*_pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(childIndex, coarserLayerIndex);
                tsl::sparse_set<unsigned int>::const_iterator it3;
                for (it3 = childIndexParentList.begin(); it3 != childIndexParentList.end(); it3++) {
                    unsigned int const childParentIndex = (*it3);

                    if (coarserPivotLayerNeighborsList.find(childParentIndex) == coarserPivotLayerNeighborsList.end()) {
                        flag_addChild = false;
                        break;
                    }
                }
                if (!flag_addChild) continue;

                // check that all parents of query are coarse neighbors to childIndex
                // if (_pivotLayers->coarseNeighborsAvailable) {
                // if (!_pivotLayers->coarseNeighborsAvailable) {
                //     printf("ERROR: COARSE NEIGHBORS NOT AVAILABLE??");
                // }

                tsl::sparse_set<unsigned int> const& childCoarseNeighborsList =
                    (*_pivotLayers)[finerLayerIndex].get_coarsePivotLayerNeighbors(childIndex);
                tsl::sparse_set<unsigned int>::const_iterator it4;
                for (it4 = parentsList.begin(); it4 != parentsList.end(); it4++) {
                    unsigned int const parentIndex = (*it4);

                    if (childCoarseNeighborsList.find(parentIndex) == childCoarseNeighborsList.end()) {
                        flag_addChild = false;
                        break;
                    }
                }
                if (!flag_addChild) continue;
                // }

                localPotentialNeighborsList.insert(childIndex);
            }
        }

// update potential neighbors list
#pragma omp critical(stage3_potentialNeighborsListUpdate)
        { potentialFineLayerPivotNeighborsList.insert(localPotentialNeighborsList.begin(), localPotentialNeighborsList.end()); }
    }

    // compute distance to all remaining fine pivots
    // if no lune -> automatic neighbor. else, add to activeChildrenList
    std::vector<unsigned int> potentialNeighborsVector(potentialFineLayerPivotNeighborsList.begin(), potentialFineLayerPivotNeighborsList.end());
#pragma omp parallel num_threads(_numThreads)
    {
        // create local sparse matrix and local neighbors list
        SearchSparseMatrixLocal localSparseMatrix = _sparseMatrix->getLocalCopy();
        tsl::sparse_set<unsigned int> localFinerPivotLayerNeighborsList;
        tsl::sparse_set<unsigned int> localActiveChildrenList;

#pragma omp for schedule(static)
        for (int it3 = 0; it3 < (int)potentialNeighborsVector.size(); it3++) {
            unsigned int const potentialNeighborIndex = potentialNeighborsVector[it3];
            float const childRadius = (*_pivotLayers)[finerLayerIndex].get_pivotRadius(potentialNeighborIndex);
            float const distance = _sparseMatrix->getDistance_LocalAndGlobal(queryIndex, potentialNeighborIndex, localSparseMatrix);

            // if no lune, add to neighbor list
            if ((_luneRadius(distance, queryRadius, childRadius) <= 0) || (_luneRadius(distance, childRadius, queryRadius) <= 0)) {
                localFinerPivotLayerNeighborsList.insert(potentialNeighborIndex);
            } else {
                localActiveChildrenList.insert(potentialNeighborIndex);
            }
        }

#pragma omp critical(stage3_updateFinerLayerNeighbors)
        { finerPivotLayerNeighborsList.insert(localFinerPivotLayerNeighborsList.begin(), localFinerPivotLayerNeighborsList.end()); }

#pragma omp critical(stage3_updateActiveChildren)
        { activeChildrenList.insert(localActiveChildrenList.begin(), localActiveChildrenList.end()); }

#pragma omp critical(stage3_updateSparseMatrix)
        { _sparseMatrix->updateSparseMatrixWithLocal(localSparseMatrix); }
    }

    // sort activeChildrenList to stay same order
    std::vector<unsigned int> activeChildrenListVectorCopy(activeChildrenList.begin(), activeChildrenList.end());
    std::sort(activeChildrenListVectorCopy.begin(), activeChildrenListVectorCopy.end());
    activeChildrenList.clear();
    activeChildrenList.insert(activeChildrenListVectorCopy.begin(), activeChildrenListVectorCopy.end());

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[3] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[3] += (dEnd_stage - dStart_stage);
    return;
}




/**
 * @brief Pivot-Exemplar Interactions: check that exemplars are coarse neighbors with all parents of Q.
 *
 * @param queryStruct
 */
void OnlineSearch::stage3(QueryStructOnlineSearch& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    // initializations
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int> const& coarserPivotLayerNeighborsList = queryStruct.pivotLayerNeighborsListMap[coarserLayerIndex];
    tsl::sparse_set<unsigned int> const& parentsList = queryStruct.ancestorIndicesMap[coarserLayerIndex];
    tsl::sparse_set<unsigned int>& activeChildrenList = queryStruct.activeChildrenList;
    tsl::sparse_set<unsigned int>& finerPivotLayerNeighborsList = queryStruct.pivotLayerNeighborsListMap[finerLayerIndex];

    //===============================================================================
    // collect all children of coarse neighbors
    //===============================================================================
    tsl::sparse_set<unsigned int> potentialFineLayerPivotNeighborsList{};

    // iterate through each coarse neighbor
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarserPivotLayerNeighborsList.begin(); it1 != coarserPivotLayerNeighborsList.end(); it1++) {
        unsigned int const coarserNeighborIndex = (*it1);
        tsl::sparse_set<unsigned int> const& coarseNeighborChildrenList =
            (*_pivotLayers)[coarserLayerIndex].get_descendantPivotIndices(coarserNeighborIndex, finerLayerIndex);

        // add all children to list
        potentialFineLayerPivotNeighborsList.insert(coarseNeighborChildrenList.begin(), coarseNeighborChildrenList.end());
    }

    //===============================================================================
    // check all parents of fine pivots are coarse neighbors of query, and vice versa
    //===============================================================================
    tsl::sparse_set<unsigned int>::const_iterator it2;
    for (it2 = potentialFineLayerPivotNeighborsList.begin(); it2 != potentialFineLayerPivotNeighborsList.end(); it2++) {
        unsigned int const childIndex = (*it2);
        bool flag_addChild = true;

        // ensure query has all parents of child as coarse neighbor
        tsl::sparse_set<unsigned int> const& childParentList =
            (*_pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(childIndex, coarserLayerIndex);
        tsl::sparse_set<unsigned int>::const_iterator it3;
        for (it3 = childParentList.begin(); it3 != childParentList.end(); it3++) {
            unsigned int const childParentIndex = (*it3);

            if (coarserPivotLayerNeighborsList.find(childParentIndex) == coarserPivotLayerNeighborsList.end()) {
                flag_addChild = false;
                break;
            }
        }
        if (!flag_addChild) continue;

        // ensure parents of query are coarse neighbors of child index
        if (_pivotLayers->coarseNeighborsAvailable) {
            tsl::sparse_set<unsigned int> const& childCoarseNeighborsList = (*_pivotLayers)[finerLayerIndex].get_coarsePivotLayerNeighbors(childIndex);
            tsl::sparse_set<unsigned int>::const_iterator it4;
            for (it4 = parentsList.begin(); it4 != parentsList.end(); it4++) {
                unsigned int const parentIndex = (*it4);

                if (childCoarseNeighborsList.find(parentIndex) == childCoarseNeighborsList.end()) {
                    flag_addChild = false;
                    break;
                }
            }
            if (!flag_addChild) continue;
        }
        

        // now, childIndex either becomes neighbor or more consideration is needed
        float const childRadius = (*_pivotLayers)[finerLayerIndex].get_pivotRadius(childIndex);
        float const distance = getDistance(queryIndex, childIndex);

        if (_luneRadius(distance, queryRadius, childRadius) <= 0 || _luneRadius(distance, childRadius, queryRadius) <= 0) {
            finerPivotLayerNeighborsList.insert(childIndex);
        } else {
            activeChildrenList.insert(childIndex);
        }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[3] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[3] += (dEnd_stage - dStart_stage);
}