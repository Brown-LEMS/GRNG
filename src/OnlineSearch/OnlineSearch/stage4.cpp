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
 * @brief Pivot-Nediated Exemplar-Exemplar Interactions: Check pivots for interference
 *
 * @param queryStruct
 */
void OnlineSearch::stage4(QueryStructOnlineSearch& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int>& activeChildrenList = queryStruct.activeChildrenList;
    RankOrderedList const& rankOrderedCoarsePivotsList = queryStruct.rankOrderedCoarsePivotsList;
    RankOrderedList& rankOrderedFinePivotsList = queryStruct.rankOrderedFinePivotsList;

    // consider each pivot in activeChildrenList as a potential neighbor
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeChildrenList.begin(); it1 != activeChildrenList.end(); /* increment in loop */) {
        unsigned int const activeChildIndex = (*it1);
        float const activeChildRadius = (*_pivotLayers)[finerLayerIndex].get_pivotRadius(activeChildIndex);
        float const distance12 = getDistance(queryIndex, activeChildIndex);
        bool flag_erase = false;

        // check for interference by coarse pivots in rank ordered list
        std::vector<IndexDistanceStruct>::const_iterator it2;
        for (it2 = rankOrderedCoarsePivotsList.list.begin(); it2 != rankOrderedCoarsePivotsList.list.end(); it2++) {
            unsigned int const interferingPivotIndex = (*it2).index;
            float const distance13 = (*it2).distance;

            if (interferingPivotIndex == activeChildIndex || interferingPivotIndex == queryIndex) continue;

            // lune side 1 check, can break since ranked list
            if (distance13 >= _luneRadius(distance12, queryRadius, activeChildRadius)) break;

            // lune side 2 check if distance is available
            std::pair<bool, float> const distance23_pair = getDistanceIfAvailable(activeChildIndex, interferingPivotIndex);
            if (distance23_pair.first) {
                if (distance23_pair.second < _luneRadius(distance12, activeChildRadius, queryRadius)) {
                    flag_erase = true;
                    break;
                }
            }
        }

        if (flag_erase) {
            it1 = activeChildrenList.erase(it1);
        } else {
            it1++;
            rankOrderedFinePivotsList.add(activeChildIndex, distance12);
        }
    }

    rankOrderedFinePivotsList.sort();
    rankOrderedFinePivotsList.resize(25);  // research here

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[4] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[4] += (dEnd_stage - dStart_stage);
}