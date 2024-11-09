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
// 02-03-2022
#include "offline-build.hpp"

#include <iostream>

#include "constant-radius-pivot-selection.hpp"
#include "memory-usage.h"

OfflineBuild::OfflineBuild(std::shared_ptr<PivotLayerHierarchy> const& pivotLayers)
    : _pivotLayers(pivotLayers),
      _sparseMatrix(pivotLayers->get_sparseMatrix()),
      _pivotRadiusVector(pivotLayers->get_pivotRadiusVector()),
      _numberOfPivotLayers(_pivotRadiusVector.size()),
      _datasetSize(_sparseMatrix->get_datasetSize()),
      _luneType(pivotLayers->get_luneType()) {
    _setLuneType(_luneType);
    _pivotLayers->coarseNeighborsAvailable = true;
    return;
}

/**
 * @brief perform the offline build on pivotLayerHierarchy
 *
 */
void OfflineBuild::build(bool verbose) {
    if (verbose) std::cout << "\nBegin Offline Index Construction..." << std::endl;
    initializeMetrics();
    std::chrono::high_resolution_clock::time_point tStart_build = std::chrono::high_resolution_clock::now();

    // pivot selection
    if (verbose) std::cout << "* Begin Constant Radius Pivot Selection..." << std::endl;
    unsigned long long int dStart_pivotSelection = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_pivotSelection = std::chrono::high_resolution_clock::now();
    ConstantRadiusPivotSelection pivotSelectionObject(_datasetSize, _pivotRadiusVector, _sparseMatrix);
    std::vector<PivotSelectionStruct> selectedPivotIndicesAndRadii = pivotSelectionObject.selectPivots();
    std::chrono::high_resolution_clock::time_point tEnd_pivotSelection = std::chrono::high_resolution_clock::now();
    unsigned long long int dEnd_pivotSelection = get_distanceComputationCount();
    _time_pivotSelection = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_pivotSelection - tStart_pivotSelection);
    _distanceCount_pivotSelection = (dEnd_pivotSelection - dStart_pivotSelection);

    // number of pivots per layer
    if (verbose) std::cout << "    + Number of Pivots Per Layer: ";
    _numberOfPivotsPerLayer.resize(_numberOfPivotLayers);
    for (int layerIndex = 0; layerIndex < _numberOfPivotLayers; layerIndex++) {
        _numberOfPivotsPerLayer[layerIndex] = selectedPivotIndicesAndRadii[layerIndex].pivotIndices->size();
        if (verbose) std::cout << selectedPivotIndicesAndRadii[layerIndex].pivotIndices->size() << ", ";
    }
    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "* End Constant Radius Pivot Selection!" << std::endl;

    // top-down incremental index construction
    if (verbose) std::cout << "* Begin Top-Down Index Construction..." << std::endl;
    for (int layerIndex = 0; layerIndex < _numberOfPivotLayers; layerIndex++) {
        if (verbose) std::cout << "    + Begin Layer " << layerIndex << " Construction..." << std::endl;
        unsigned long long int dStart_layer = get_distanceComputationCount();
        std::chrono::high_resolution_clock::time_point tStart_layer = std::chrono::high_resolution_clock::now();

        tsl::sparse_set<unsigned int>& pivotIndices = *selectedPivotIndicesAndRadii[layerIndex].pivotIndices;
        tsl::sparse_map<unsigned int, float>& pivotRadii = *selectedPivotIndicesAndRadii[layerIndex].pivotRadii;

        // brute force construction
        if (layerIndex == 0) {
            (*_pivotLayers).setPivotIndicesAndRadii(0, pivotIndices, pivotRadii);
            (*_pivotLayers)[0].computeBruteForcePivotLayerGraph();

        } else {
            tsl::sparse_set<unsigned int>::const_iterator it1;
            for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
                unsigned int const queryIndex = (*it1);
                float const queryRadius = pivotRadii[queryIndex];
                QueryStructOfflineBuild queryStruct(layerIndex, queryIndex, queryRadius);

                // check if query is already a pivot in the above layer
                queryStruct.flag_isCoarsePivot = false;
                tsl::sparse_set<unsigned int> const& currentPivotIndices = *(*_pivotLayers)[layerIndex].get_pivotIndices_ptr();
                if (currentPivotIndices.find(queryIndex) != currentPivotIndices.end()) {
                    queryStruct.flag_isCoarsePivot = true;
                }

                // find neighbors of query in a top-down approach
                for (int coarserLayerIndex = 0; coarserLayerIndex < layerIndex; coarserLayerIndex++) {
                    queryStruct.newCoarserLayer(coarserLayerIndex);

                    stage0(queryStruct);  // Preprocessing

                    // needed to find parents in each layer
                    // but skipping to coarser layer, working off work previously done
                    if (queryStruct.flag_isCoarsePivot) {
                        if (coarserLayerIndex < layerIndex - 1) continue;
                    }

                    stage1(queryStruct);  // Pivot-Pivot Interactions

                    stage2(queryStruct);  // Query-Pivot Interactions

                    stage3_multithreading(queryStruct);  // Pivot-Exemplar Interactions

                    stage4(queryStruct);  // Pivot-Mediated Exemplar-Exemplar Interactions

                    stage5(queryStruct);  // Exemplar-Mediated Exemplar-Exemplar Interactions

                    stage6_multithreading(queryStruct);  // Verification of Query Links
                }

                stage7_multithreading(queryStruct);  // Verification of Exisiting Links

                updateDataStructures(queryStruct);
            }
        }

        std::chrono::high_resolution_clock::time_point tEnd_layer = std::chrono::high_resolution_clock::now();
        unsigned long long int dEnd_layer = get_distanceComputationCount();
        _time_layers[layerIndex] = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_layer - tStart_layer);
        _distanceCount_layers[layerIndex] = (dEnd_layer - dStart_layer);
        if (verbose) std::cout << "    + End Layer " << layerIndex << " Construction!" << std::endl;
    }
    if (verbose) std::cout << "* End Top-Down Index Construction!" << std::endl;

    std::chrono::high_resolution_clock::time_point tEnd_build = std::chrono::high_resolution_clock::now();
    unsigned long long int dEnd_build = get_distanceComputationCount();
    _time_build = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_build - tStart_build);
    _distanceCount_build = dEnd_build;
    if (verbose) std::cout << "End Offline Index Construction! \n" << std::endl;
    if (verbose) printMetrics();
}

void OfflineBuild::selectPivots() {
    printf("\nBegin Pivot Selection...");
    ConstantRadiusPivotSelection pivotSelectionObject(_datasetSize, _pivotRadiusVector, _sparseMatrix);
    std::vector<PivotSelectionStruct> selectedPivotIndicesAndRadii = pivotSelectionObject.selectPivots();

    // number of pivots per layer
    printf("    + Number of Pivots Per Layer: \n");
    _numberOfPivotsPerLayer.resize(_numberOfPivotLayers);
    for (int layerIndex = 0; layerIndex < _numberOfPivotLayers; layerIndex++) {
        _numberOfPivotsPerLayer[layerIndex] = selectedPivotIndicesAndRadii[layerIndex].pivotIndices->size();
        std::cout << selectedPivotIndicesAndRadii[layerIndex].pivotIndices->size() << ", ";
    }
    std::cout << std::endl;
    std::cout << "* End Constant Radius Pivot Selection!" << std::endl;

    return;    
}


//===================================================================
// UPDATE DATA STRUCTURES
//  - note -> order matters:
//  - addPivotToLayer to inialize all data structures
//  - update parents/children (impacts umax and dmax) NEXT
//  - remove links (smaller set to recompute umax from)
//  - add links
//  - add coarse pivots
//===================================================================
void OfflineBuild::updateDataStructures(QueryStructOfflineBuild& queryStruct) {
    unsigned long long int const dStart_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_stage = std::chrono::high_resolution_clock::now();

    // initializations
    unsigned int const queryIndex = queryStruct.queryIndex;
    float const queryRadius = queryStruct.queryRadius;
    int const pivotLayerIndex = queryStruct.pivotLayerIndex;
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap = queryStruct.ancestorIndicesMap;
    LinksMap const& invalidatedPivotLinks = queryStruct.invalidatedPivotLinks;
    tsl::sparse_set<unsigned int> const& finerLayerPivotNeighborList = queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex];
    tsl::sparse_set<unsigned int> const& coarseNeighborList = queryStruct.pivotLayerNeighborsListMap[pivotLayerIndex - 1];

    // add pivot to layer and update structures
    (*_pivotLayers).addPivotToLayer(pivotLayerIndex, queryIndex, queryRadius);

    // update ancestor/descendant relations
    (*_pivotLayers).updateAncestors(pivotLayerIndex, queryIndex, ancestorIndicesMap);

    // remove the links invalidated by query
    (*_pivotLayers).removePivotLayerNeighbors(pivotLayerIndex, invalidatedPivotLinks.map);

    // add links to query
    (*_pivotLayers).addPivotLayerNeighbors(pivotLayerIndex, queryIndex, finerLayerPivotNeighborList);

    // add the coarse layer neighbors of query
    (*_pivotLayers).addCoarsePivotLayerNeighbors(pivotLayerIndex, queryIndex, coarseNeighborList);

    unsigned long long int const dEnd_stage = get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tEnd_stage = std::chrono::high_resolution_clock::now();
    _time_stages[8] = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_stage - tStart_stage);
    _distanceCount_stages[8] = (dEnd_stage - dStart_stage);
    return;
}

void OfflineBuild::printMetrics() {
    std::cout << "===================================================" << std::endl;
    std::cout << "===== Offline Index Construction Statistics  ======" << std::endl;
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
    std::cout << "Number of Threads: " << _numThreads << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "  Timings: (s)" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Total Build (s): " << _time_build.count() << std::endl;
    std::cout << "Pivot Selection (s): " << _time_pivotSelection.count() << std::endl;
    std::cout << "Layers: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _time_layers[i].count() << ", ";
    }
    std::cout << std::endl;
    std::cout << "Stages: ";
    for (int i = 0; i < 9; i++) {
        std::cout << _time_stages[i].count() << ", ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "  Distance Computations: " << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Total Build: " << _distanceCount_build << std::endl;
    std::cout << "Pivot Selection: " << _distanceCount_pivotSelection << std::endl;
    std::cout << "Layers: ";
    for (int i = 0; i < _numberOfPivotLayers; i++) {
        std::cout << _distanceCount_layers[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Stages: ";
    for (int i = 0; i < 9; i++) {
        std::cout << _distanceCount_stages[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "  Memory Consumption: " << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Current RSS: " << std::endl;
    std::cout << "    (B): " << (float)(getCurrentRSS()) << std::endl;
    std::cout << "    (KB): " << (float)(getCurrentRSS() / (double)(1 << 10)) << std::endl;
    std::cout << "    (MB): " << (float)(getCurrentRSS() / (double)(1 << 20)) << std::endl;
    std::cout << "    (GB): " << (float)(getCurrentRSS() / (double)(1 << 30)) << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "===================================================" << std::endl;
}