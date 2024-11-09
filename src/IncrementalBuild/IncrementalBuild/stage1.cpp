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

/**
 * @brief Pivot-Pivot Interactions. Collect pivots that are neighbors to all parents of Q.
 *
 * @param queryStruct
 */
void IncrementalBuild::stage1(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    tsl::sparse_set<unsigned int> const& parentIndices = queryStruct.ancestorIndicesMap[coarserLayerIndex];
    tsl::sparse_set<unsigned int> const& coarserLayerNeighborsList =
        queryStruct.pivotLayerNeighborsListMap[coarserLayerIndex];
    tsl::sparse_set<unsigned int>& activeCoarsePivotsList = queryStruct.activeCoarsePivotsList;

    // if neighbors available from previous iteration
    if (!coarserLayerNeighborsList.empty()) {
        activeCoarsePivotsList = coarserLayerNeighborsList;
        
    } else {

        //===============================================================================
        // collect potential coarser neighbors
        //===============================================================================
        if (queryStruct.flag_isCoarsePivot) { /// OUT OF ORDER.
            // if query is a pivot in the coarse layer
            // then we can use the neighbors from above layer
            activeCoarsePivotsList = (*_pivotLayers)[coarserLayerIndex].get_pivotLayerNeighbors(queryIndex);
            activeCoarsePivotsList.insert(queryIndex);

        } else {
            // collect neighbors of parents
            tsl::sparse_set<unsigned int>::const_iterator it1;
            for (it1 = parentIndices.begin(); it1 != parentIndices.end(); it1++) {
                unsigned int const parentIndex = (*it1);
                tsl::sparse_set<unsigned int> const& parentNeighborsList =
                    (*_pivotLayers)[coarserLayerIndex].get_pivotLayerNeighbors(parentIndex);
                activeCoarsePivotsList.insert(parentIndex);
                activeCoarsePivotsList.insert(parentNeighborsList.begin(), parentNeighborsList.end());
            }
        }

        //===============================================================================
        // remove pivots not neighbors to all parents
        //===============================================================================
        tsl::sparse_set<unsigned int>::iterator it2;
        for (it2 = activeCoarsePivotsList.begin(); it2 != activeCoarsePivotsList.end(); /* increment in loop */) {
            unsigned int const activeCoarsePivotIndex = (*it2);
            bool flag_erase = false;

            // not removing parents as they become neighbors
            if (parentIndices.find(activeCoarsePivotIndex) != parentIndices.end()) {
                it2++;
                continue;
            }

            // check that activeCoarsePivotIndex is neighbor to all parents
            tsl::sparse_set<unsigned int> const& activeCoarsePivotNeighborsList =
                (*_pivotLayers)[coarserLayerIndex].get_pivotLayerNeighbors(activeCoarsePivotIndex);
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = parentIndices.begin(); it3 != parentIndices.end(); it3++) {
                unsigned int const parentIndex = (*it3);

                if (activeCoarsePivotNeighborsList.find(parentIndex) == activeCoarsePivotNeighborsList.end()) {
                    flag_erase = true;
                }
            }

            // remove from list if applicable
            if (flag_erase) {
                it2 = activeCoarsePivotsList.erase(it2);
            } else {
                it2++;
            }
        }
    }

    // create ranked list from coarse pivots
    RankOrderedList& rankOrderedCoarsePivotsList = queryStruct.rankOrderedCoarsePivotsList;
    RankOrderedList& rankOrderedErodedLuneCoarsePivotsList = queryStruct.rankOrderedErodedLuneCoarsePivotsList;

    tsl::sparse_set<unsigned int>::const_iterator it4;
    for (it4 = activeCoarsePivotsList.begin(); it4 != activeCoarsePivotsList.end(); it4++) {
        unsigned int const activeCoarsePivotIndex = (*it4);
        float const activeCoarsePivotRadius =
            (*_pivotLayers)[coarserLayerIndex].get_pivotRadius(activeCoarsePivotIndex);

        float const distance = getDistance(queryIndex, activeCoarsePivotIndex);
        float const erodedDistance = _luneRadius(distance, queryRadius, activeCoarsePivotRadius);

        // add to ranked lists
        rankOrderedCoarsePivotsList.add(activeCoarsePivotIndex, distance);
        rankOrderedErodedLuneCoarsePivotsList.add(activeCoarsePivotIndex, erodedDistance);
    }

    rankOrderedCoarsePivotsList.sort();
    rankOrderedErodedLuneCoarsePivotsList.sort();

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[1] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[1] += (dEnd_stage - dStart_stage);

    return;
}