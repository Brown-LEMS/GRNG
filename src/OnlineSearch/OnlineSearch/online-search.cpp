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
 * @brief Initialize Online Search object by providing constructed hierarchy
 * 
 * @param pivotLayers 
 */
OnlineSearch::OnlineSearch(std::shared_ptr<PivotLayerHierarchy> const& pivotLayers)
    : _pivotLayers(pivotLayers),
      _datasetSize(pivotLayers->get_datasetSize()),
      _numberOfPivotLayers(pivotLayers->get_numberOfPivotLayers()),
      _luneType(pivotLayers->get_luneType()) {
    _setLuneType(_luneType);
}

/**
 * @brief perform online search over a test set
 * 
 * @param testPointer 
 * @param testsetSize 
 * @param dimension 
 * @param verbose 
 * @param searchNeighbors returns RNG neighbors of each query
 */
void OnlineSearch::search(float* const& testPointer, unsigned int const testsetSize, unsigned int const dimension, bool verbose,
                          tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& searchNeighbors) {
    if (verbose) printf("Begin Online Search for %u queries...\n", testsetSize);
    initializeMetrics();
    _testsetSize = testsetSize;

    if (verbose) printf("Coarse Neighbors Access?: %u\n", _pivotLayers->coarseNeighborsAvailable);

    // initializing sparse matrix for search
    if (dimension != _pivotLayers->get_dimension()) {
        printf("Testset dimension does not match dataset. Aborting search!\n");
        return;
    }
    _sparseMatrix = std::make_shared<SearchSparseMatrix>(testPointer, testsetSize, dimension, _pivotLayers->get_sparseMatrix());

    // find neighbors in bottom layer
    int const pivotLayerIndex = _numberOfPivotLayers - 1;
    float const queryRadius = 0.0f;

    // find neighbors of each query
    searchNeighbors.clear();

    unsigned long long int dStart_search = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_search = std::chrono::high_resolution_clock::now();
    for (int testIndex = 0; testIndex < testsetSize; testIndex++) {
        searchNeighbors[testIndex] = tsl::sparse_set<unsigned int>{};
        unsigned int const queryIndex = _datasetSize + testIndex;

        if (verbose) printf("queryIndex: %u, neighbors: ", testIndex);
        _search(pivotLayerIndex, queryIndex, queryRadius, searchNeighbors[testIndex]);
        if (verbose) printSet(searchNeighbors[testIndex]);
    }
    std::chrono::high_resolution_clock::time_point tEnd_search = std::chrono::high_resolution_clock::now();
    unsigned long long int dEnd_search = get_distanceComputationCount();
    _time_averageQuery = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_search - tStart_search) / testsetSize;
    _distanceCount_averageQuery = ((float)(dEnd_search - dStart_search)) / testsetSize;
    if (verbose) printMetrics();
    return;
}

/**
 * @brief find the BG Neighbors of (queryIndex,queryRadius) in layer pivotLayerIndex
 * 
 * @param pivotLayerIndex 
 * @param queryIndex 
 * @param queryRadius 
 * @param queryNeighbors BG neighbors updated here
 */
void OnlineSearch::_search(int const pivotLayerIndex, unsigned int const queryIndex, float const queryRadius,
                           tsl::sparse_set<unsigned int>& queryNeighbors) {
    QueryStructOnlineSearch queryStruct(pivotLayerIndex, queryIndex, queryRadius);

    if (pivotLayerIndex == 0) {  // single-layer build by brute force
        findNeighborsBruteForce(queryStruct);

    } else {  // multi-layer in hierarchical fashion
        for (int coarserLayerIndex = 0; coarserLayerIndex < pivotLayerIndex; coarserLayerIndex++) {
            queryStruct.newCoarserLayer(coarserLayerIndex);

            stage0(queryStruct);

            stage1(queryStruct);

            stage2(queryStruct);

            stage3_multithreading(queryStruct);

            stage4(queryStruct);

            stage5(queryStruct);

            stage6_multithreading(queryStruct);
        }
    }

    queryNeighbors = queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex];
    return;
}

void OnlineSearch::printMetrics() {
    printf("===================================================\n");
    printf("===          Online Search Statistics          ====\n");
    printf("===================================================\n");
    printf("\n");
    printf("Number of Queries: %u\n", _testsetSize);
    printf("Number of Threads: %u\n", _numThreads);
    printf("\n");
    printf("---------------------------------------------------\n");
    printf("  Timings: (s)\n");
    printf("---------------------------------------------------\n");
    printf("\n");
    printf("Average Query (ms): %.3f\n", _time_averageQuery.count() * 1000);
    printf("\n");
    printf("---------------------------------------------------\n");
    printf("  Distance Computations: \n");
    printf("---------------------------------------------------\n");
    printf("\n");
    printf("Average Query: %.2f\n", _distanceCount_averageQuery);
    printf("\n");
    printf("Time Stages: \n");
    for (int i = 0; i < 7; i++) {
        printf("%.6f,", 1000*_time_stages[i].count());
    }
    printf("\n");
    printf("Distance Stages: \n");    
    for (int i = 0; i < 7; i++) {
        printf("%llu,", _distanceCount_stages[i]);
    }
    printf("\n");

    printf("---------------------------------------------------\n");
    printf("\n");
    printf("===================================================\n");
}

void OnlineSearch::printSet(tsl::sparse_set<unsigned int> const& set) {
    printf("{");
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}
