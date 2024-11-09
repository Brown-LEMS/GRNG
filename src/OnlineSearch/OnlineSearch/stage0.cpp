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
#include "online-search.hpp"

/**
 * @brief Preprocessing Stage. Find all parents of Q within the layer.
 *
 * @param queryStruct
 */
void OnlineSearch::stage0(QueryStructOnlineSearch& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    tsl::sparse_set<unsigned int>& parentsList = queryStruct.ancestorIndicesMap[coarserLayerIndex];
    RankOrderedList rankOrderedTopLayerPivotsList;

    //===============================================================================
    // create the list of potential parents
    //===============================================================================
    tsl::sparse_set<unsigned int> activeParentList;
    if (coarserLayerIndex == 0) {
        // all top layer pivots
        tsl::sparse_set<unsigned int> const& topLayerPivotIndices = *(*_pivotLayers)[coarserLayerIndex].get_pivotIndices_ptr();
        activeParentList.insert(topLayerPivotIndices.begin(), topLayerPivotIndices.end());

    } else {
        // collect all grandparents
        tsl::sparse_set<unsigned int> const& grandParentsList = queryStruct.ancestorIndicesMap[coarserLayerIndex - 1];

        // collect all children of grandparents (potential parents)
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = grandParentsList.begin(); it1 != grandParentsList.end(); it1++) {
            unsigned int const grandParentIndex = (*it1);
            tsl::sparse_set<unsigned int> const& grandParentPivotChildrenList =
                (*_pivotLayers)[coarserLayerIndex - 1].get_descendantPivotIndices(grandParentIndex, coarserLayerIndex);

            activeParentList.insert(grandParentPivotChildrenList.begin(), grandParentPivotChildrenList.end());
        }
    }

    //===============================================================================
    // find the true parents from the potential
    //===============================================================================
    tsl::sparse_set<unsigned int>::const_iterator it2;
    for (it2 = activeParentList.begin(); it2 != activeParentList.end(); it2++) {
        unsigned int const activeParentIndex = (*it2);
        float const activeParentPivotRadius = (*_pivotLayers)[coarserLayerIndex].get_pivotRadius(activeParentIndex);
        float const distance = getDistance(queryIndex, activeParentIndex);

        // parent if within radius
        if (distance <= (activeParentPivotRadius - queryRadius)) {
            parentsList.insert(activeParentIndex);
        }

        // for use in stage6
        if (coarserLayerIndex == 0) {
            rankOrderedTopLayerPivotsList.add(activeParentIndex, distance);
        }
    }

    // updating data structures
    if (coarserLayerIndex == 0) {
        rankOrderedTopLayerPivotsList.sort();
        queryStruct.rankOrderedTopLayerPivotsList = rankOrderedTopLayerPivotsList;
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[0] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[0] += (dEnd_stage - dStart_stage);
}
