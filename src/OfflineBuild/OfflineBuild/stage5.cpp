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
// 2022-02-04
#include "offline-build.hpp"

/**
 * @brief Exemplar-mediated Exemplar-Exemplar Interactions: check exemplar for interference
 *
 * @param queryStruct
 */
void OfflineBuild::stage5(QueryStructOfflineBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const coarserLayerIndex = queryStruct.coarserLayerIndex;
    int const finerLayerIndex = coarserLayerIndex + 1;
    tsl::sparse_set<unsigned int>& activeChildrenList = queryStruct.activeChildrenList;
    RankOrderedList const& rankOrderedCoarsePivotsList = queryStruct.rankOrderedCoarsePivotsList;
    RankOrderedList& rankOrderedFinePivotsList = queryStruct.rankOrderedFinePivotsList;

    unsigned int fanoutThreshold = 25;

    // consider each pivot in activeChildrenList as a potential neighbor
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = activeChildrenList.begin(); it1 != activeChildrenList.end(); /* increment in loop */) {
        unsigned int const activeChildIndex = (*it1);
        float const activeChildRadius = (*_pivotLayers)[finerLayerIndex].get_pivotRadius(activeChildIndex);
        float const distance12 = getDistance(queryIndex, activeChildIndex);
        bool flag_erase = false;

        // check for interference by fine-layer pivots
        std::vector<IndexDistanceStruct>::const_iterator it2;
        for (it2 = rankOrderedFinePivotsList.list.begin(); it2 != rankOrderedCoarsePivotsList.list.end(); it2++) {
            unsigned int const interferingPivotIndex = (*it2).index;
            float const distance13 = (*it2).distance;

            // lune side1 check
            // ranked list, so power to break
            if (distance13 >= _luneRadius(distance12, queryRadius, activeChildRadius)) break;

            // lune side2 check
            std::pair<bool, float> const distance23_pair = getDistanceIfAvailable(activeChildIndex, interferingPivotIndex);
            if (distance23_pair.first) {
                if (distance23_pair.second < _luneRadius(distance12, activeChildRadius, queryRadius)) {
                    flag_erase = true;
                    break;
                }
            }
        }

        // test neighbors of activeChildIndex for interference
        if (!flag_erase) {
            unsigned int fanoutCounter = 0;

            tsl::sparse_set<unsigned int> const& activeChildNeighborsList =
                (*_pivotLayers)[finerLayerIndex].get_pivotLayerNeighbors(activeChildIndex);

            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = activeChildNeighborsList.begin(); it3 != activeChildNeighborsList.end() &&
                                                         fanoutCounter < fanoutThreshold;
                 it3++) {
                unsigned int const interferingPivotIndex = (*it3);

                std::pair<bool, float> const distance13_pair = getDistanceIfAvailable(queryIndex, interferingPivotIndex);
                std::pair<bool, float> const distance23_pair = getDistanceIfAvailable(activeChildIndex, interferingPivotIndex);

                if (distance13_pair.first && distance23_pair.first) {
                    if (distance13_pair.second < _luneRadius(distance12, queryRadius, activeChildRadius)) {
                        if (distance23_pair.second < _luneRadius(distance12, activeChildRadius, queryRadius)) {
                            flag_erase = true;
                            break;
                        }
                    }
                }
                fanoutCounter++;
            }
        }

        // erase or increment potential neighbor
        if (flag_erase) {
            it1 = activeChildrenList.erase(it1);
        } else {
            it1++;
        }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[5] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[5] += (dEnd_stage - dStart_stage);
}