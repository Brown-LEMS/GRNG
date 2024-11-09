#ifndef QueryStructIncrementalBuild_hpp
#define QueryStructIncrementalBuild_hpp
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

#include <tsl/sparse_map.h>
#include <tsl/sparse_set.h>

#include "links-map.hpp"
#include "rank-ordered-list.hpp"

/**
 * @brief Struct to hold all temporary datastructures for individual query during OfflineBuild
 *
 */
struct QueryStructIncrementalBuild {
  public:
    QueryStructIncrementalBuild(){};
    QueryStructIncrementalBuild(unsigned int const queryIndex) : queryIndex(queryIndex){};
    ~QueryStructIncrementalBuild(){};

    // QUERY INFORMATION
    unsigned int const queryIndex = 0;
    int pivotLayerIndex = 0;
    float queryRadius = 0.0f;

    /**
     * @brief Adding query to Layer layerIndex. Updating and resetting all necessary datastructures
     *
     * @param layerIndex
     */
    inline void newPivotLayer(int const layerIndex) {
        pivotLayerIndex = layerIndex;
        queryRadius = queryRadiusMap[pivotLayerIndex];
        ancestorIndicesMap = ancestorIndicesMapMap[pivotLayerIndex];
        descendantIndicesMap = descendantIndicesMapMap[pivotLayerIndex];
        pivotLayerNeighborsListMap = tsl::sparse_map<int, tsl::sparse_set<unsigned int>> {};
        for (int i = 0; i <= pivotLayerIndex; i++) {
            pivotLayerNeighborsListMap[i] = tsl::sparse_set<unsigned int> {};
        }
        finerCoarseNeighborsList = tsl::sparse_set<unsigned int> {};
        activeCoarsePivotsList = tsl::sparse_set<unsigned int>{};
        rankOrderedCoarsePivotsList = RankOrderedList();
        rankOrderedErodedLuneCoarsePivotsList = RankOrderedList();
        activeChildrenList = tsl::sparse_set<unsigned int>{};
        rankOrderedFinePivotsList = RankOrderedList();
        invalidatedPivotLinks = LinksMap();
        invalidatedFineCoarseNeighborsLinks = LinksMap();
    }

    // layer information
    int coarserLayerIndex = 0;
    inline void newCoarserLayer(int const layer) {
        coarserLayerIndex = layer;
        activeCoarsePivotsList.clear();
        rankOrderedCoarsePivotsList.clear();
        rankOrderedErodedLuneCoarsePivotsList.clear();
        rankOrderedFinePivotsList.clear();
        activeChildrenList.clear();
    }

    // layer where query first becomes a pivot
    int orphanLayerID = -1;
    tsl::sparse_map<int, float> queryRadiusMap{};
    int flag_isCoarsePivot = false;

    // maps of ancestors of query in each layer its a pivot
    tsl::sparse_map<int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> ancestorIndicesMapMap{};
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> ancestorIndicesMap{};

    // maps of descendants of query in each layer its a pivot
    tsl::sparse_map<int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> descendantIndicesMapMap{};
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> descendantIndicesMap{};

    // maps of neighbors of query for each layer its a pivot
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> pivotLayerNeighborsListMap{};

    // ADDED FOR VIRTUAL/COARSE NEIGHBORS 2022-06-15
    tsl::sparse_set<unsigned int> finerCoarseNeighborsList{};

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
    LinksMap invalidatedFineCoarseNeighborsLinks = LinksMap();

    // indices of importance to query
    bool topLayerPivotsAdded = false;
    tsl::sparse_set<unsigned int> _queryImportantIndices{};
    void updateImportantIndices() {

        // add ancestors
        tsl::sparse_map<int,tsl::sparse_set<unsigned int>>::const_iterator it1;
        for (it1 = ancestorIndicesMap.begin(); it1 != ancestorIndicesMap.end(); it1++) {
            tsl::sparse_set<unsigned int> const& ancestors = (*it1).second;
            _queryImportantIndices.insert(ancestors.begin(),ancestors.end());
        }

        // add descendants
        tsl::sparse_map<int,tsl::sparse_set<unsigned int>>::const_iterator it2;
        for (it2 = descendantIndicesMap.begin(); it2 != descendantIndicesMap.end(); it2++) {
            tsl::sparse_set<unsigned int> const& descenants = (*it2).second;
            _queryImportantIndices.insert(descenants.begin(),descenants.end());
        }

        // add neighbors in current layer (ADD COARSER LAYER)
        _queryImportantIndices.insert(pivotLayerNeighborsListMap[pivotLayerIndex].begin(),pivotLayerNeighborsListMap[pivotLayerIndex].end());
        if (pivotLayerIndex > 0) {
            _queryImportantIndices.insert(pivotLayerNeighborsListMap[pivotLayerIndex-1].begin(),pivotLayerNeighborsListMap[pivotLayerIndex-1].end());
        }

        // added finer layer coarse neighbors**
        _queryImportantIndices.insert(finerCoarseNeighborsList.begin(), finerCoarseNeighborsList.end());

        // storing all pairwise top-layer distances
        if (pivotLayerIndex == 0) {
            std::vector<IndexDistanceStruct>::const_iterator it3;
            for (it3 = rankOrderedTopLayerPivotsList.list.begin(); it3 != rankOrderedTopLayerPivotsList.list.end(); it3++) {
                _queryImportantIndices.insert((*it3).index);
            }
        }

        // add all top layer pivots
        // if (!topLayerPivotsAdded) {
        //     std::vector<IndexDistanceStruct>::const_iterator it3;
        //     for (it3 = rankOrderedTopLayerPivotsList.list.begin(); it3 != rankOrderedTopLayerPivotsList.list.end(); it3++) {
        //         _queryImportantIndices.insert((*it3).index);
        //     }
        //     topLayerPivotsAdded = true;
        // }

    }
};

#endif  // QueryStructIncrementalBuild_hpp