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
// 2022-05-29
#include "incremental-build.hpp"

#include <iostream>

#include "memory-usage.h"

IncrementalBuild::IncrementalBuild(std::shared_ptr<PivotLayerHierarchy> const& pivotLayers)
    : _pivotLayers(pivotLayers),
      _sparseMatrix(pivotLayers->get_sparseMatrix()),
      _pivotRadiusVector(pivotLayers->get_pivotRadiusVector()),
      _numberOfPivotLayers(_pivotRadiusVector.size()),
      _datasetSize(_sparseMatrix->get_datasetSize()),
      _luneType(pivotLayers->get_luneType()) {
    _setLuneType(_luneType);
    _cacheAll = _sparseMatrix->getCacheAllStatus();
    _pivotLayers->coarseNeighborsAvailable = true;

}

/**
 * @brief perform incremental build using
 *
 * @param verbose
 */
void IncrementalBuild::build(bool verbose) {
    if (verbose) std::cout << "\nBegin Incremental Index Construction..." << std::endl;
    initializeMetrics();
    std::chrono::high_resolution_clock::time_point tStart_build = std::chrono::high_resolution_clock::now();
    if (verbose) printf("%-15s%-20s%-20s\n", "queryIndex", "averageDistances", "averageTime(ms)");

    // incremental build
    int statsBatchSize = 100;
    unsigned long long int dStart_query = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_query = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < _datasetSize; queryIndex++) {
        // add query to hierarchy
        incrementalUpdate(queryIndex);

        // update average query statistics
        if ((queryIndex + 1) % statsBatchSize == 0) {
            std::chrono::high_resolution_clock::time_point tEnd_query = std::chrono::high_resolution_clock::now();
            unsigned long long int dEnd_query = get_distanceComputationCount();
            _time_averageQuery = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_query - tStart_query) / statsBatchSize;
            _distanceCount_averageQuery = ((float)(dEnd_query - dStart_query)) / statsBatchSize;

            if (verbose) {
                printf("%-20u%-20.2f%-20.3f\n", queryIndex + 1, _distanceCount_averageQuery, _time_averageQuery.count() * 1000);
            }
            dStart_query = get_distanceComputationCount();
            tStart_query = std::chrono::high_resolution_clock::now();
        }
    }

    std::chrono::high_resolution_clock::time_point tEnd_build = std::chrono::high_resolution_clock::now();
    unsigned long long int dEnd_build = get_distanceComputationCount();
    _time_build = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_build - tStart_build);
    _distanceCount_build = dEnd_build;

    _numberOfPivotsPerLayer.clear();
    _numberOfPivotsPerLayer = (*_pivotLayers).get_numberOfPivotsPerLayer();
    (*_pivotLayers).getDegreeStats(_numberLinksPerLayer,_averageDegreePerLayer,_minDegreePerLayer,_maxDegreePerLayer);
    _currentRSS = (float)(getCurrentRSS());


    if (verbose) std::cout << "End Incremental Index Construction! \n" << std::endl;
    if (verbose) printMetrics();
}

void IncrementalBuild::incrementalUpdate(unsigned int queryIndex) {

    // initialize structures for new query
    if (!_cacheAll) {
        _sparseMatrix->setNewQuery(queryIndex); // sparse matrix to temp cache all query distances
    }
    QueryStructIncrementalBuild queryStruct(queryIndex);

    // find layers where query is to be added
    // find all ancestors for all layers above each layer query is added
    // find all descendants for all layers below each layer query is added
    stage0(queryStruct);

    // add query to each layer
    for (int pivotLayerIndex = queryStruct.orphanLayerID; pivotLayerIndex < _numberOfPivotLayers; pivotLayerIndex++) {
        // std::cout << "  * pivotLayerIndex: " << pivotLayerIndex << std::endl;
        queryStruct.newPivotLayer(pivotLayerIndex);

        // adding pivot by brute force to the layer
        if (pivotLayerIndex == 0) {
            findNeighborsBruteForce(queryStruct);

            findInterferenceBruteForce(queryStruct);

        } else {
            // // did we already add query to the coarser layer?
            // queryStruct.flag_isCoarsePivot = false;
            // if (pivotLayerIndex > queryStruct.orphanLayerID) {
            //     queryStruct.flag_isCoarsePivot = true;
            // }

            // find neighbors of (Q,r_l) in layer l
            // working top-down
            for (int coarserLayerIndex = 0; coarserLayerIndex < pivotLayerIndex; coarserLayerIndex++) {
                queryStruct.newCoarserLayer(coarserLayerIndex);

                // if flag_coarserLayerPivot, skip to layer l-1 and use neighbors of (Q,r_{l-1}) as potential coarse neighbors
                // if (queryStruct.flag_isCoarsePivot) {
                //     if (coarserLayerIndex < pivotLayerIndex - 1) continue;
                // }

                // collect potential coarse neighbors as neighbors. else get neighbors
                stage1(queryStruct);

                // find coarse neighbors
                stage2(queryStruct);

                // collect children of coarse neighbors as potential neighbors
                stage3_multithreading(queryStruct);

                //
                stage4(queryStruct);

                stage5(queryStruct);

                // exhaustively check links for interference
                stage6_multithreading(queryStruct);
            }

            // check if Q invalidates any links
            stage7_multithreading(queryStruct);
        }

        // if coarse pivot, Q might be a coarse neighbor of fine pivots. or it might interfere with current ones
        if (pivotLayerIndex < _numberOfPivotLayers - 1) {

            // Stage8: find fine pivots that Q is a coarse neighbor of 
            stage8_multithreading(queryStruct);

            // Stage9: find fine pivot- coarse pivot links that Q invalidates
            stage9_multithreading(queryStruct);
        }

        updateDataStructures(queryStruct);
    }

    // add query distances to sparse matrix 
    if (!_cacheAll) {
        _sparseMatrix->updateQueryMatrix(queryStruct._queryImportantIndices);
    }

    return;
}

void IncrementalBuild::updateDataStructures(QueryStructIncrementalBuild& queryStruct) {
    unsigned long long int const dStart_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();

    //  initializations
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap = queryStruct.ancestorIndicesMap;
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& descendantIndicesMap = queryStruct.descendantIndicesMap;
    LinksMap const& invalidatedPivotLinks = queryStruct.invalidatedPivotLinks;
    tsl::sparse_set<unsigned int> const& pivotLayerNeighborsList = queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex];


    // add pivot to layer and update structures
    (*_pivotLayers).addPivotToLayer(pivotLayerIndex, queryIndex, queryRadius);

    // update ancestor/descendant relations
    (*_pivotLayers).updateAncestors(pivotLayerIndex, queryIndex, ancestorIndicesMap); // this first for a reason. umax/vmax

    // update ancestor/descendant relations
    (*_pivotLayers).updateDescendants(pivotLayerIndex, queryIndex, descendantIndicesMap);

    // remove the links invalidated by query
    (*_pivotLayers).removePivotLayerNeighbors(pivotLayerIndex, invalidatedPivotLinks.map);

    // add links to query
    (*_pivotLayers).addPivotLayerNeighbors(pivotLayerIndex, queryIndex, pivotLayerNeighborsList);

    // remove coarse neighbors first
    if (pivotLayerIndex < _numberOfPivotLayers - 1) {
        
        // first remove coarse neighbbors. updates cMaxDescendants and all
        (*_pivotLayers).removeFinerCoarsePivotLayerNeighbors(pivotLayerIndex,queryStruct.invalidatedFineCoarseNeighborsLinks.map);

        // now add coarse neighbors
        (*_pivotLayers).addFinerCoarsePivotLayerNeighbors(pivotLayerIndex,queryIndex,queryStruct.finerCoarseNeighborsList);
    }

    if (pivotLayerIndex > 0) {
        // add the coarse layer neighbors of query
        (*_pivotLayers).addCoarsePivotLayerNeighbors(pivotLayerIndex, queryIndex, queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex-1]);
    }

    // update query list of relationship indices
    if (!_cacheAll) {
        queryStruct.updateImportantIndices();
    }

    // add AFTER query is added as top layer pivot
    if (pivotLayerIndex == 0) {
        queryStruct.rankOrderedTopLayerPivotsList.add(queryIndex,0.0f);
        queryStruct.rankOrderedTopLayerPivotsList.sort();
    }

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[10] += std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[10] += (dEnd_stage - dStart_stage);
    return;
}

void IncrementalBuild::printMetrics() {
    std::cout << "===================================================" << std::endl;
    std::cout << "=== Incremental Index Construction Statistics  ====" << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Lune Type: " << _luneType << std::endl;
    std::cout << "Number of Layers: " << _numberOfPivotLayers << std::endl;
    std::cout << "Radii Per Layer: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _pivotRadiusVector[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Pivots Per Layer: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _numberOfPivotsPerLayer[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Links Per Layer: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _numberLinksPerLayer[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Ave Degree Per Layer: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _averageDegreePerLayer[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Min Degree Per Layer: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _minDegreePerLayer[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Max Degree Per Layer: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _maxDegreePerLayer[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Number of Threads: " << _numThreads << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "  Timings: (s)" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Total Build (s): " << _time_build.count() << std::endl;
    std::cout << "Layers: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _time_layers[i].count() << ", ";
    }
    std::cout << std::endl;
    std::cout << "Stages: ";
    for (int i = 0; i < 11; i++) {
        std::cout << _time_stages[i].count() << ", ";
    }
    std::cout << std::endl;
    std::cout << "Average Query (ms): " << _time_averageQuery.count() * 1000 << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "  Distance Computations: " << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Total Build: " << _distanceCount_build << std::endl;
    std::cout << "Layers: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _distanceCount_layers[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Stages: ";
    for (int i = 0; i < 11; i++) {
        std::cout << _distanceCount_stages[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Average Query: " << _distanceCount_averageQuery << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "  Memory Consumption: " << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Current RSS: " << std::endl;
    std::cout << "    (B): " << _currentRSS << std::endl;
    std::cout << "    (KB): " << (float)(_currentRSS / (double)(1 << 10)) << std::endl;
    std::cout << "    (MB): " << (float)(_currentRSS / (double)(1 << 20)) << std::endl;
    std::cout << "    (GB): " << (float)(_currentRSS / (double)(1 << 30)) << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    printf("L, R, P, Index Distances, Index Time (s), Memory (GB), Ave. Degree, Min Degree, Max Degree\n");
    printf("%u,",_numberOfPivotLayers);
    printf("{;");
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        printf("%.6f;", _pivotRadiusVector[i]);
    } printf("},");
    printf("{;");
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        printf("%u;", _numberOfPivotsPerLayer[i]);
    } printf("},");
    printf("%llu,",_distanceCount_build);
    printf("%.6f,",_time_build.count());
    printf("%.6f,",(float)(_currentRSS / (double)(1 << 30)));


    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "===================================================" << std::endl;
}