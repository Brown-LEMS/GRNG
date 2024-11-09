#include "constant-radius-pivot-selection.hpp"

#include <iostream>

#include "constant-radius-pivot-selection-struct.hpp"

ConstantRadiusPivotSelection::ConstantRadiusPivotSelection(unsigned int const datasetSize,
                                                           std::vector<float> const& pivotRadiusVector, 
                                                           std::shared_ptr<SparseMatrix> const& sparseMatrix) : _datasetSize(datasetSize), _pivotRadiusVector(pivotRadiusVector), _numberOfPivotLayers(pivotRadiusVector.size()), _sparseMatrix(sparseMatrix) {};

//===================================================================
// Pivot Selection for Offline Build
//===================================================================
std::vector<PivotSelectionStruct> ConstantRadiusPivotSelection::selectPivots() {
    // intializing each parent layer with ConstantRadiusPivotSelectionStruct
    std::vector<ConstantRadiusPivotSelectionStruct> pivotSelectionVector;
    for (int layerIndex = 0; layerIndex < _numberOfPivotLayers - 1; layerIndex++) {
        pivotSelectionVector.push_back(ConstantRadiusPivotSelectionStruct(layerIndex,
                                                                          _pivotRadiusVector[layerIndex], _numberOfPivotLayers));
    }

    // add bottom Layer, not necesarrily RNG
    float bottomLayerRadius = _pivotRadiusVector[_numberOfPivotLayers - 1];

    // GREEDILY, INCREMENTALLY ADD PIVOTS TO LAYERS
    for (unsigned int queryIndex = 0; queryIndex < _datasetSize; queryIndex++) {
        // initializations for each pivot
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>> parentIndices;
        tsl::sparse_set<int> layerIndicesToAdd;
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>> memberIndicesToAdd;
        for (int layerIndex = 0; layerIndex < (int)pivotSelectionVector.size(); layerIndex++) {
            parentIndices[layerIndex] = tsl::sparse_set<unsigned int>{};
            memberIndicesToAdd[layerIndex] = tsl::sparse_set<unsigned int>{};
        }

        // FOR EACH PIVOT, FIND POTENTIAL PARENTS, AKA QUERY WITHIN RADIUS
        for (int layerIndex = 0; layerIndex < (int)pivotSelectionVector.size(); layerIndex++) {
            // GATHER CANDIDATE PARENT INDICES
            tsl::sparse_set<unsigned int> candidateParentIndices;
            if (layerIndex == 0) {
                candidateParentIndices = *pivotSelectionVector[layerIndex].pivotIndices;
            } else {
                // potential parents as children of parents in previous layer
                tsl::sparse_set<unsigned int> const& previousLayerParentIndices = parentIndices[layerIndex - 1];
                tsl::sparse_set<unsigned int>::const_iterator it1;
                for (it1 = previousLayerParentIndices.begin(); it1 != previousLayerParentIndices.end(); it1++) {
                    unsigned int const parentIndex = (*it1);
                    tsl::sparse_set<unsigned int> const& childIndices =
                        (*pivotSelectionVector[layerIndex - 1].pivotMemberIndices)[parentIndex];
                    candidateParentIndices.insert(childIndices.begin(), childIndices.end());
                }
            }

            // CHECK IF QUERY WITHIN RADIUS
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = candidateParentIndices.begin(); it2 != candidateParentIndices.end(); it2++) {
                unsigned int const candidateParentIndex = (*it2);
                float distance = getDistance(queryIndex, candidateParentIndex);

                if (distance <= _pivotRadiusVector[layerIndex] - bottomLayerRadius) {  // CHANGE FOR NO RNG LAYER
                    parentIndices[layerIndex].insert(candidateParentIndex);
                }
            }
        }

        // BOTTOM-UP APPROACH TO FIND QUERY AS TRUE MEMBER OF PIVOTS
        // aka, if query has any true parents, r_p - r_q
        // determines if query becomes a pivot
        int currentLayerIndex = pivotSelectionVector.size() - 1;

        // if no parents at bottom(ish) layer, becomes a pivot
        // if parents, not a pivot in any layer (but rng) so ignore loop
        if (parentIndices[currentLayerIndex].empty()) {
            layerIndicesToAdd.insert(currentLayerIndex);
            currentLayerIndex--;

            // test for true parents in above layers
            bool flag_isPivot = true;
            while (flag_isPivot && currentLayerIndex >= 0) {
                // if no potential parents, becomes a pivot, next layer
                if (parentIndices[currentLayerIndex].empty()) {
                    layerIndicesToAdd.insert(currentLayerIndex);
                    memberIndicesToAdd[currentLayerIndex].insert(queryIndex);
                } else {
                    // for potential parents, check if any are true parents
                    float parentRadius = _pivotRadiusVector[currentLayerIndex];
                    float queryRadius = _pivotRadiusVector[currentLayerIndex + 1];

                    tsl::sparse_set<unsigned int> const& candidateParentIndices = parentIndices[currentLayerIndex];
                    tsl::sparse_set<unsigned int>::const_iterator it3;
                    for (it3 = candidateParentIndices.begin(); it3 != candidateParentIndices.end(); it3++) {
                        unsigned int const candidateParentIndex = (*it3);
                        float distance = getDistance(queryIndex, candidateParentIndex);

                        if (distance <= parentRadius - queryRadius) {
                            flag_isPivot = false;
                            memberIndicesToAdd[currentLayerIndex].insert(candidateParentIndex);
                        }
                    }

                    // if no parents, query is a pivot, next layer
                    if (flag_isPivot == true) {
                        layerIndicesToAdd.insert(currentLayerIndex);
                        memberIndicesToAdd[currentLayerIndex].insert(queryIndex);
                    }
                }

                // increment to higher layer
                currentLayerIndex--;
            }

            // update query information in pivotSelectionVector
            for (int layerIndex = 0; layerIndex < (int)pivotSelectionVector.size(); layerIndex++) {
                // add pivots to layer
                if (layerIndicesToAdd.find(layerIndex) != layerIndicesToAdd.end()) {  // in
                    pivotSelectionVector[layerIndex].addPivot(queryIndex);
                }

                // add pivot members to layer
                if (layerIndex < (int)pivotSelectionVector.size() - 1) {
                    tsl::sparse_set<unsigned int>::const_iterator it4;
                    for (it4 = memberIndicesToAdd[layerIndex].begin(); it4 != memberIndicesToAdd[layerIndex].end(); it4++) {
                        unsigned int const currentMemberIndex = (*it4);
                        pivotSelectionVector[layerIndex].addPivotMember(currentMemberIndex, queryIndex);
                    }
                }
            }
        }
    }

    // create PivotSelectionStruct() for each layer, holding pivotIndices and radii
    std::vector<PivotSelectionStruct> selectedPivotIndicesAndRadii;
    for (int layerIndex = 0; layerIndex < (int)pivotSelectionVector.size(); layerIndex++) {
        std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices_ptr =
            std::make_shared<tsl::sparse_set<unsigned int>>();
        std::shared_ptr<tsl::sparse_map<unsigned int, float>> pivotRadii_ptr =
            std::make_shared<tsl::sparse_map<unsigned int, float>>();

        tsl::sparse_set<unsigned int> const& layer_PivotIndices = *pivotSelectionVector[layerIndex].pivotIndices;
        tsl::sparse_set<unsigned int>::const_iterator it5;
        for (it5 = layer_PivotIndices.begin(); it5 != layer_PivotIndices.end(); it5++) {
            unsigned int const pivotIndex = (*it5);
            pivotIndices_ptr->insert(pivotIndex);
            pivotRadii_ptr->emplace(pivotIndex, _pivotRadiusVector[layerIndex]);
        }
        selectedPivotIndicesAndRadii.push_back(PivotSelectionStruct(pivotIndices_ptr, pivotRadii_ptr));
    }

    // add bottom layer
    int layerIndex = _numberOfPivotLayers - 1;
    std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices_ptr =
        std::make_shared<tsl::sparse_set<unsigned int>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, float>> pivotRadii_ptr =
        std::make_shared<tsl::sparse_map<unsigned int, float>>();
    for (unsigned int queryIndex = 0; queryIndex < _datasetSize; queryIndex++) {
        pivotIndices_ptr->insert(queryIndex);
        pivotRadii_ptr->emplace(queryIndex, _pivotRadiusVector[layerIndex]);
    }
    selectedPivotIndicesAndRadii.push_back(PivotSelectionStruct(pivotIndices_ptr, pivotRadii_ptr));

    return selectedPivotIndicesAndRadii;
}