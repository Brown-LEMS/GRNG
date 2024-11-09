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
 * @brief get neighbors of query within the layer. used for finding neighbors on top layer and debugging.
 * 
 * @param queryStruct 
 */
void OnlineSearch::findNeighborsBruteForce(QueryStructOnlineSearch& queryStruct) {
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    tsl::sparse_set<unsigned int>& pivotLayerNeighbors = queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex];

    // consider each pivot as potential neighbor of queryIndex
    tsl::sparse_set<unsigned int> const& pivotIndices = *((*_pivotLayers)[pivotLayerIndex].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1, it2;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const potentialNeighborIndex = (*it1);
        float const potentialNeighborRadius = (*_pivotLayers)[pivotLayerIndex].get_pivotRadius(potentialNeighborIndex);
        float const distance12 = getDistance(queryIndex, potentialNeighborIndex);
        float const queryLinkDistance = _luneRadius(distance12, queryRadius, potentialNeighborRadius);
        float const potentialNeighborLinkDistance = _luneRadius(distance12, potentialNeighborRadius, queryRadius);
        bool flag_isNeighbor = true;

        // no lune?
        if (queryLinkDistance <= 0 || potentialNeighborLinkDistance <= 0) {
            pivotLayerNeighbors.insert(potentialNeighborIndex);
            continue;
        }

        // consider each pivot as interfering with link(queryIndex,potentialNeighborIndex)
        for (it2 = pivotIndices.begin(); it2 != pivotIndices.end(); it2++) {
            unsigned int const interferingPivotIndex = (*it2);
            if (interferingPivotIndex == potentialNeighborIndex || interferingPivotIndex == queryIndex) continue;

            // check for lune interference
            float const distance23 = getDistance(potentialNeighborIndex, interferingPivotIndex);
            if (distance23 < potentialNeighborLinkDistance) {
                float const distance13 = getDistance(queryIndex, interferingPivotIndex);
                if (distance13 < queryLinkDistance) {
                    flag_isNeighbor = false;
                    break;
                }
            }
        }

        if (flag_isNeighbor) {
            pivotLayerNeighbors.insert(potentialNeighborIndex);
        }
    }

    return;
}