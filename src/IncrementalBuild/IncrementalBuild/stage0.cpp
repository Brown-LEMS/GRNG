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

#include <iostream>

#include "incremental-build.hpp"
#include "rank-ordered-list.hpp"  // temp

// want to find each layer that query belongs to, and parents in each layer
void IncrementalBuild::stage0(QueryStructIncrementalBuild& queryStruct) {
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();
    unsigned long long int const dStart_stage = get_distanceComputationCount();

    unsigned int const queryIndex = queryStruct.queryIndex;
    tsl::sparse_map<int, float>& queryRadiusMap = queryStruct.queryRadiusMap;
    int& orphanLayerID = queryStruct.orphanLayerID;
    orphanLayerID = _numberOfPivotLayers - 1;
    tsl::sparse_map<int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>& ancestorIndicesMap = queryStruct.ancestorIndicesMapMap;
    tsl::sparse_map<int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>& descendantIndicesMap = queryStruct.descendantIndicesMapMap;
    RankOrderedList& rankOrderedTopLayerPivotsList = queryStruct.rankOrderedTopLayerPivotsList;
    rankOrderedTopLayerPivotsList.clear();

    // working bottom-up. adding query to layerIndex=0..L-1
    ancestorIndicesMap.clear();
    bool flag_finished = false;
    for (int layerIndex = _numberOfPivotLayers - 1; layerIndex >= 0 && !flag_finished; layerIndex--) {
        float const queryRadius = _pivotRadiusVector[layerIndex];
        queryRadiusMap[layerIndex] = queryRadius;

        if (layerIndex == 0) continue;  // no parents to be found
        ancestorIndicesMap[layerIndex] = tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{};

        // find parents of query in coarserLayerIndex = 0..layerIndex - 1
        // if query in layerIndex has parents in layerIndex - 1, query is not added to layerIndex - 1
        for (int coarserLayerIndex = 0; coarserLayerIndex < layerIndex; coarserLayerIndex++) {
            ancestorIndicesMap[layerIndex][coarserLayerIndex] = tsl::sparse_set<unsigned int>{};

            // if query is an orphan, won't have parents in this layer
            if (coarserLayerIndex >= orphanLayerID) {
                ancestorIndicesMap[layerIndex][coarserLayerIndex].insert(queryIndex);
                continue;  //
            }

            // collect list of potential parents
            tsl::sparse_set<unsigned int> potentialParentsList{};
            if (layerIndex == _numberOfPivotLayers - 1) {  // first time through

                if (coarserLayerIndex == 0) {  // add all top layer pivots
                    potentialParentsList = *((*_pivotLayers)[coarserLayerIndex].get_pivotIndices_ptr());

                } else {  // add all children of grandparents
                    tsl::sparse_set<unsigned int> const& grandParentsList = ancestorIndicesMap[layerIndex][coarserLayerIndex - 1];
                    tsl::sparse_set<unsigned int>::const_iterator it1;
                    for (it1 = grandParentsList.begin(); it1 != grandParentsList.end(); it1++) {
                        unsigned int const grandParent = (*it1);
                        tsl::sparse_set<unsigned int> const& grandParentChildren =
                            (*_pivotLayers)[coarserLayerIndex - 1].get_descendantPivotIndices(grandParent, coarserLayerIndex);
                        potentialParentsList.insert(grandParentChildren.begin(), grandParentChildren.end());
                    }
                }

            } else {  // already know parents of query in this layer when query has smaller radius.

                // parents of (Q,r_Q1) is a subset of parents of (Q,r_Q2) for r_Q1 > r_Q2
                potentialParentsList = ancestorIndicesMap[layerIndex + 1][coarserLayerIndex];
            }

            // find true parents
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = potentialParentsList.begin(); it2 != potentialParentsList.end(); it2++) {
                unsigned int const potentialParent = (*it2);
                float const potentialParentRadius = (*_pivotLayers)[coarserLayerIndex].get_pivotRadius(potentialParent);
                float const distance = _sparseMatrix->getDistance(queryIndex, potentialParent);

                // true parent condition
                if (distance <= (potentialParentRadius - queryRadius)) {
                    ancestorIndicesMap[layerIndex][coarserLayerIndex].insert(potentialParent);
                }

                // add to top layer ranked list if first time through
                if ((layerIndex == _numberOfPivotLayers - 1) && (coarserLayerIndex == 0)) {
                    rankOrderedTopLayerPivotsList.add(potentialParent, distance);
                }
            }

            // if no parents, then query will be a pivot in that layer!
            if (ancestorIndicesMap[layerIndex][coarserLayerIndex].empty()) {
                orphanLayerID = coarserLayerIndex;

                // consider all lower layer pivots having queryIndex as a parent
                for (int previousPivotLayerIndex = layerIndex; previousPivotLayerIndex < _numberOfPivotLayers; previousPivotLayerIndex++) {
                    // query becomes pivot in all layers below too
                    for (int previousCoarserLayerIndex = orphanLayerID; previousCoarserLayerIndex < previousPivotLayerIndex;
                         previousCoarserLayerIndex++) {
                        ancestorIndicesMap[previousPivotLayerIndex][previousCoarserLayerIndex].insert(queryIndex);
                    }
                }

                // if have parents in layerIndex-1, then query does not become pivot in that layer. exit!
            } else {
                if (coarserLayerIndex == layerIndex - 1) {
                    if (!ancestorIndicesMap[layerIndex][coarserLayerIndex].empty()) {
                        flag_finished = true;
                        break;
                    }
                }
            }
        }
    }

    // sort the rank ordered list
    // if (orphanLayerID == 0) {
    //     rankOrderedTopLayerPivotsList.add(queryIndex, 0);
    // }
    rankOrderedTopLayerPivotsList.sort();

    // now, collect existing children of queryIndex for each layer its being added to
    // this is a top-down situation
    //
    // collecting "intersecting" pivot indices
    // i.e. pivot pi s.t. d(Q,pi) <= ri + Q. aka pivot domains intersect
    //
    // start with highest layer, l = orphanLayer, largest radius for Q.
    //      find intersecting indices from lbar = 0..L-1
    //      check all top layer pivots or children of above intersecting
    //      if lbar > l, then check for descendants
    // go down to l = orphanLayer+1 .. L-2
    //      find intersecting indices from lbar = l..L-1
    //      check intersecting from query with larger radius (subset. pivot must be closer since radius smaller)
    //      if lbar > l, then check for descendants

    // collect intersecting pivot indices for each layer above orphan
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> intersectingPivotIndices;
    descendantIndicesMap.clear();

    // for each "coarse" layer that query becomes a pivot
    for (int pivotLayerIndex = orphanLayerID; pivotLayerIndex < _numberOfPivotLayers - 1; pivotLayerIndex++) {
        float const queryRadius = queryRadiusMap[pivotLayerIndex];
        descendantIndicesMap[pivotLayerIndex] = tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{};

        if (pivotLayerIndex == orphanLayerID) {  // only on first go around
            for (int layerIndex = 0; layerIndex < _numberOfPivotLayers; layerIndex++) {
                if (layerIndex > pivotLayerIndex) {
                    descendantIndicesMap[pivotLayerIndex][layerIndex] = tsl::sparse_set<unsigned int>{};
                }

                intersectingPivotIndices[layerIndex] = tsl::sparse_set<unsigned int>{};
                if (layerIndex == 0) {
                    // all top layer pivots as potentially intersecting
                    tsl::sparse_set<unsigned int> const& topLayerPivots = *(*_pivotLayers)[0].get_pivotIndices_ptr();
                    tsl::sparse_set<unsigned int>::const_iterator it1;
                    for (it1 = topLayerPivots.begin(); it1 != topLayerPivots.end(); it1++) {
                        unsigned int const topLayerPivot = (*it1);
                        float const topLayerPivotRadius = (*_pivotLayers)[layerIndex].get_pivotRadius(topLayerPivot);
                        float const distance = _sparseMatrix->getDistance(queryIndex, topLayerPivot);

                        // unfortunately cannot break since pivot radii not factored in
                        if (distance < (topLayerPivotRadius + queryRadius)) {
                            intersectingPivotIndices[layerIndex].insert(topLayerPivot);
                        }
                    }

                } else {
                    // children of intersecting as possibly intersecting
                    tsl::sparse_set<unsigned int> potentialIntersectingPivotsList{};
                    tsl::sparse_set<unsigned int> const& intersectingParentsList = intersectingPivotIndices[layerIndex - 1];
                    tsl::sparse_set<unsigned int>::const_iterator it1;
                    for (it1 = intersectingParentsList.begin(); it1 != intersectingParentsList.end(); it1++) {
                        unsigned int const parentIndex = (*it1);
                        tsl::sparse_set<unsigned int> const& potentialIntersectingChildrenList =
                            (*_pivotLayers)[layerIndex - 1].get_descendantPivotIndices(parentIndex, layerIndex);
                        potentialIntersectingPivotsList.insert(potentialIntersectingChildrenList.begin(), potentialIntersectingChildrenList.end());
                    }

                    // find intersecting from potentiall intersecting
                    tsl::sparse_set<unsigned int>::const_iterator it2;
                    for (it2 = potentialIntersectingPivotsList.begin(); it2 != potentialIntersectingPivotsList.end(); it2++) {
                        unsigned int const potentialIntersectingIndex = (*it2);
                        float const potentialIntersectingRadius = (*_pivotLayers)[layerIndex].get_pivotRadius(potentialIntersectingIndex);
                        float const distance = getDistance(queryIndex, potentialIntersectingIndex);

                        // unfortunately cannot break since pivot radii not factored in
                        if (distance <= (potentialIntersectingRadius + queryRadius)) {
                            intersectingPivotIndices[layerIndex].insert(potentialIntersectingIndex);

                            // is this intersecting pivot a descendant?
                            if (layerIndex > pivotLayerIndex) {
                                if (distance <= (queryRadius - potentialIntersectingRadius)) {
                                    descendantIndicesMap[pivotLayerIndex][layerIndex].insert(potentialIntersectingIndex);
                                }
                            }
                        }
                    }
                }
            }

        } else {
            // intersecting indices as a subset of previous intersecting indices with larger r_Q
            for (int layerIndex = pivotLayerIndex + 1; layerIndex < _numberOfPivotLayers; layerIndex++) {
                descendantIndicesMap[pivotLayerIndex][layerIndex] = tsl::sparse_set<unsigned int>{};

                // get and clear previous intersecting indices
                tsl::sparse_set<unsigned int> const potentialIntersecting = intersectingPivotIndices[layerIndex];
                intersectingPivotIndices[layerIndex].clear();

                // find new intersecting from previous intersecting
                tsl::sparse_set<unsigned int>::const_iterator it1;
                for (it1 = potentialIntersecting.begin(); it1 != potentialIntersecting.end(); it1++) {
                    unsigned int const potentialIntersectingIndex = (*it1);
                    float const potentialIntersectingRadius = (*_pivotLayers)[layerIndex].get_pivotRadius(potentialIntersectingIndex);
                    float const distance = getDistance(queryIndex, potentialIntersectingIndex);

                    // unfortunately cannot break since pivot radii not factored in
                    if (distance <= (potentialIntersectingRadius + queryRadius)) {
                        intersectingPivotIndices[layerIndex].insert(potentialIntersectingIndex);

                        // is this intersecting pivot a descendant?
                        if (layerIndex > pivotLayerIndex) {
                            if (distance <= queryRadius - potentialIntersectingRadius) {
                                descendantIndicesMap[pivotLayerIndex][layerIndex].insert(potentialIntersectingIndex);
                            }
                        }
                    }
                }
            }
        }
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[0] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[0] += (dEnd_stage - dStart_stage);

    return;
}
