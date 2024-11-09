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
 * @brief Query-Pivot Interactions. Collect pivots that are Coarse neighbors of Q.
 *
 * @param queryStruct
 */
void OnlineSearch::stage2(QueryStructOnlineSearch& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    tsl::sparse_set<unsigned int>& coarserPivotLayerNeighborsList = queryStruct.pivotLayerNeighborsListMap[coarserLayerIndex];
    RankOrderedList& rankOrderedCoarsePivotsList = queryStruct.rankOrderedCoarsePivotsList;
    RankOrderedList& rankOrderedErodedLuneCoarsePivotsList = queryStruct.rankOrderedErodedLuneCoarsePivotsList;

    if (!coarserPivotLayerNeighborsList.empty()) return;  // already have neighbors

    // iterate through rankOrderedErodedLuneCoarsePivotsList as potential neighbors
    std::vector<IndexDistanceStruct>::const_iterator it1;
    for (it1 = rankOrderedErodedLuneCoarsePivotsList.list.begin(); it1 != rankOrderedErodedLuneCoarsePivotsList.list.end(); it1++) {
        unsigned int const activeCoarsePivotIndex = (*it1).index;
        float const activeCoarsePivotRadius = (*_pivotLayers)[coarserLayerIndex].get_pivotRadius(activeCoarsePivotIndex);
        float const distance12 = getDistance(queryIndex, activeCoarsePivotIndex);
        float const distance12_eroded = (*it1).distance;
        bool flag_isNeighbor = true;

        // no lune, no interference possible
        if (distance12_eroded <= 0) {
            coarserPivotLayerNeighborsList.insert(activeCoarsePivotIndex);
            continue;
        }

        // iterate through rankOrderedCoarsePivotsList as potentially interfering
        std::vector<IndexDistanceStruct>::const_iterator it2;
        for (it2 = rankOrderedCoarsePivotsList.list.begin(); it2 != rankOrderedCoarsePivotsList.list.end(); it2++) {
            unsigned int const interferingPivotIndex = (*it2).index;
            float const distance13 = (*it2).distance;

            // first side lune interference check.
            // ranked list, can break the search
            if (distance13 >= distance12_eroded) break;

            // second side lune interference check
            float const distance23 = getDistance(activeCoarsePivotIndex, interferingPivotIndex);
            if (distance23 < _luneRadius(distance12, activeCoarsePivotRadius, queryRadius)) {
                flag_isNeighbor = false;
                break;
            }
        }

        if (flag_isNeighbor) {
            coarserPivotLayerNeighborsList.insert(activeCoarsePivotIndex);
        }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[2] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[2] += (dEnd_stage - dStart_stage);
}