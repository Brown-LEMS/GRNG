#ifndef QueryStructOfflineBuild_hpp
#define QueryStructOfflineBuild_hpp
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

#include <tsl/sparse_map.h>
#include <tsl/sparse_set.h>

#include "links-map.hpp"
#include "rank-ordered-list.hpp"

/**
 * @brief Struct to hold all temporary datastructures for individual query during OfflineBuild
 *
 */
struct QueryStructOfflineBuild {
  public:
    // CONSTRUCTORS/DESTRUCTORS
    QueryStructOfflineBuild(){};
    QueryStructOfflineBuild(int const pivotLayerIndex, unsigned int const queryIndex, float const queryRadius) 
    : pivotLayerIndex(pivotLayerIndex), queryIndex(queryIndex), queryRadius(queryRadius) {
        for (int layerIndex = 0; layerIndex < pivotLayerIndex; layerIndex++) {
            ancestorIndicesMap[layerIndex] = {};
        }

        for (int layerIndex = 0; layerIndex <= pivotLayerIndex; layerIndex++) {
            pivotLayerNeighborsListMap[layerIndex] = {};
        }
    };
    ~QueryStructOfflineBuild(){};

    // QUERY INFORMATION
    int const pivotLayerIndex = 0;
    unsigned int const queryIndex = 0;
    float const queryRadius = 0.0f;

    // layer information
    int coarserLayerIndex = 0;
    inline void newCoarserLayer(int const layer) {
        coarserLayerIndex = layer;
        activeCoarsePivotsList.clear();
        rankOrderedCoarsePivotsList.clear();
        rankOrderedErodedLuneCoarsePivotsList.clear();
        rankOrderedFinePivotsList.clear();
        activeChildrenList.clear();
        invalidatedPivotLinks.clear();
    }

    //================================================
    // DATA STRUCTURES
    //================================================

    bool flag_isCoarsePivot = false;

    // maps for each layer
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> ancestorIndicesMap =
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{};
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> pivotLayerNeighborsListMap =
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{};

    // rank ordered list to all top layer pivots
    RankOrderedList rankOrderedTopLayerPivotsList = RankOrderedList();

    // ranked list of distances to coarser layer pivots (from stage 1)
    tsl::sparse_set<unsigned int> activeCoarsePivotsList = tsl::sparse_set<unsigned int>{};
    RankOrderedList rankOrderedCoarsePivotsList = RankOrderedList();
    RankOrderedList rankOrderedErodedLuneCoarsePivotsList = RankOrderedList();

    // fine layer potential neighbors and potential interfering
    tsl::sparse_set<unsigned int> activeChildrenList = tsl::sparse_set<unsigned int>{};
    RankOrderedList rankOrderedFinePivotsList = RankOrderedList();

    // links invalidated by the addition of the query
    LinksMap invalidatedPivotLinks = LinksMap();
};

#endif  // QueryStructOfflineBuild_hpp