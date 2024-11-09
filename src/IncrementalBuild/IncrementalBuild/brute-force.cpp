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
 * @brief find neighbors of queryIndex by brute force
 * 
 * @param queryStruct 
 */
void IncrementalBuild::findNeighborsBruteForce(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

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


    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[0] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[0] += (dEnd_stage - dStart_stage);
    return;
}

/**
 * @brief Find links invalidated by queryIndex by brute force
 * 
 * @param queryStruct 
 */
void IncrementalBuild::findInterferenceBruteForce(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    unsigned int const queryIndex = queryStruct.queryIndex;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    LinksMap& invalidatedPivotLinks = queryStruct.invalidatedPivotLinks;
    LinksMap checkedLinks;

    // consider the links of each pivot
    tsl::sparse_set<unsigned int> const& pivotIndices = *((*_pivotLayers)[pivotLayerIndex].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1, it2;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex1 = (*it1);
        float const pivotRadius1 = (*_pivotLayers)[pivotLayerIndex].get_pivotRadius(pivotIndex1);
        tsl::sparse_set<unsigned int> const& pivotIndex1Neighbors = (*_pivotLayers)[pivotLayerIndex].get_pivotLayerNeighbors(pivotIndex1);

        // consider query as interfering with link(pivotIndex1,pivotIndex2)
        for (it2 = pivotIndex1Neighbors.begin(); it2 != pivotIndex1Neighbors.end(); it2++) {
            unsigned int const pivotIndex2 = (*it2);
            float const pivotRadius2 = (*_pivotLayers)[pivotLayerIndex].get_pivotRadius(pivotIndex2);
            float const distance12 = getDistance(pivotIndex1, pivotIndex2);
            bool flag_interference = false;

            // check for lune interference
            float const distance13 = getDistance(pivotIndex1, queryIndex);
            if (distance13 < _luneRadius(distance12, pivotRadius1, pivotRadius2)) {
                float const distance23 = getDistance(pivotIndex2, queryIndex);
                if (distance23 < _luneRadius(distance12, pivotRadius2, pivotRadius1)) {
                    flag_interference = true;
                }
            }

            // query inhabits lune(index1,index2), so remove link
            if (flag_interference) {
                invalidatedPivotLinks.checkAndAddLink(pivotIndex1,pivotIndex2);
            }
        }
    }


    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[0] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[0] += (dEnd_stage - dStart_stage);
    return;
}