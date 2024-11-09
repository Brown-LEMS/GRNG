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
// 01-27-2022

#include "pivot-layer-hierarchy.hpp"

#include <iostream>

/**
 * @brief Construct a new PivotLayerHierarchy, initializes SparseMatrix
 *
 * @param pivotRadiusVector
 * @param luneType
 * @param verbose
 * @param dataPointer
 * @param datasetSize
 * @param dimension
 */
PivotLayerHierarchy::PivotLayerHierarchy(std::vector<float> pivotRadiusVector, std::string luneType, float* dataPointer,
                                         unsigned int const datasetSize, unsigned int const dimension) {
    _pivotRadiusVector = pivotRadiusVector;
    _numberOfPivotLayers = pivotRadiusVector.size();
    _luneType = luneType;
    _setLuneType(luneType);

    _sparseMatrix = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension);
    initialize_hierarchy();
}

PivotLayerHierarchy::PivotLayerHierarchy(std::vector<float> const pivotRadiusVector, std::string luneType, float* dataPointer,
                                         unsigned int const datasetSize, unsigned int const dimension, bool cacheAll) {
    _pivotRadiusVector = pivotRadiusVector;
    _numberOfPivotLayers = pivotRadiusVector.size();
    _luneType = luneType;
    _setLuneType(luneType);

    _sparseMatrix = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension, cacheAll);
    initialize_hierarchy();
}

/**
 * @brief Construct a new PivotLayerHierarchy given SparseMatrix
 *
 * @param pivotRadiusVector
 * @param luneType
 * @param sparseMatrix
 */
PivotLayerHierarchy::PivotLayerHierarchy(std::vector<float> const pivotRadiusVector, std::string luneType, std::shared_ptr<SparseMatrix> sparseMatrix)

{
    _pivotRadiusVector = pivotRadiusVector;
    _numberOfPivotLayers = _pivotRadiusVector.size();
    _luneType = luneType;

    _sparseMatrix = sparseMatrix;
    initialize_hierarchy();
}

/**
 * @brief construct via existing pivotLayerHiearchy, usually loaded via serialization
 *
 * @param ref
 */
PivotLayerHierarchy::PivotLayerHierarchy(PivotLayerHierarchy const& ref) {
    _pivotRadiusVector = ref.get_pivotRadiusVector();
    _numberOfPivotLayers = _pivotRadiusVector.size();
    _luneType = ref.get_luneType();

    // copy the sparse matrix
    _sparseMatrix = ref.get_sparseMatrix();

    // copy each pivot layer
    _pivotLayers = ref.get_pivotLayers();

    // set lune type (defaults to Berk)
    if (_luneType == "Cole") {
        _luneIndex = 1;
    }
}

/**
 * @brief initialize vector of pivot layers, need sparseMatrix first
 *
 */
void PivotLayerHierarchy::initialize_hierarchy() {
    _pivotLayers.clear();
    for (int layerIndex = 0; layerIndex < _numberOfPivotLayers; layerIndex++) {
        _pivotLayers.push_back(std::make_shared<PivotLayer>(layerIndex, _numberOfPivotLayers, _sparseMatrix, _luneType, _verbose));
    }
}

/**
 * @brief Add all pivots to the layer
 *
 * @param layerIndex
 * @param pivotIndices
 * @param pivotRadii
 */
void PivotLayerHierarchy::setPivotIndicesAndRadii(int const layerIndex, tsl::sparse_set<unsigned int> const& pivotIndices,
                                                  tsl::sparse_map<unsigned int, float> const& pivotRadii) {
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        float const pivotRadius = pivotRadii.at(pivotIndex);
        (*_pivotLayers[layerIndex]).addPivotToLayer(pivotIndex, pivotRadius);
        (*_pivotLayers[layerIndex]).initializePivotInLayer(pivotIndex);
    }
}

/**
 * @brief Add a single pivot to the layer
 *
 * @param layerIndex
 * @param pivotIndex
 * @param pivotRadius
 */
void PivotLayerHierarchy::addPivotToLayer(int const layerIndex, unsigned int const pivotIndex, float const pivotRadius) {
    (*_pivotLayers[layerIndex]).addPivotToLayer(pivotIndex, pivotRadius);
    (*_pivotLayers[layerIndex]).initializePivotInLayer(pivotIndex);
}

/**
 * @brief Add all given neighbors of pivotIndex1 to the layer
 *
 * @param layerIndex
 * @param pivotIndex1
 * @param pivotLayerNeighbors
 */
void PivotLayerHierarchy::addPivotLayerNeighbors(int const layerIndex, unsigned int const pivotIndex1,
                                                 tsl::sparse_set<unsigned int> const& pivotLayerNeighbors) {
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotLayerNeighbors.begin(); it1 != pivotLayerNeighbors.end(); it1++) {
        unsigned int const pivotIndex2 = (*it1);
        if (pivotIndex2 == pivotIndex1) continue;

        (*_pivotLayers[layerIndex]).addLinkToLayer(pivotIndex1, pivotIndex2);

        if (layerIndex > 0) {
            updateUmaxDescendantMaxLinkDistance_(layerIndex, pivotIndex1);
            updateUmaxDescendantMaxLinkDistance_(layerIndex, pivotIndex2);
        }
    }
}

/**
 * @brief Remove all given links from the layer
 *
 * @param layerIndex
 * @param invalidatedPivotLinks
 */
void PivotLayerHierarchy::removePivotLayerNeighbors(int const layerIndex,
                                                    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& invalidatedPivotLinks) {
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = invalidatedPivotLinks.begin(); it1 != invalidatedPivotLinks.end(); it1++) {
        unsigned int const pivotIndex1 = (*it1).first;
        tsl::sparse_set<unsigned int> const& neighbors = (*it1).second;

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = neighbors.begin(); it2 != neighbors.end(); it2++) {
            unsigned int const pivotIndex2 = (*it2);

            // remove link from the layer
            (*_pivotLayers[layerIndex]).removeLinkFromLayer(pivotIndex1, pivotIndex2);

            if (layerIndex > 0) {
                // compute umax descendants if necessary
                checkAncestorsToRecomputeUmaxDescendant(layerIndex, pivotIndex1, pivotIndex2);
                checkAncestorsToRecomputeUmaxDescendant(layerIndex, pivotIndex2, pivotIndex1);

                // compute umax max of descendants if necessary
                checkToRecomputeUmaxMaximumOfDescendant(layerIndex, pivotIndex1, pivotIndex2);
            }
        }
    }
}

/**
 * @brief add all ancestors of pivotIndex, and add it as descendant of all ancestors
 *
 * @param layerIndex
 * @param pivotIndex
 * @param ancestorIndicesMap
 */
void PivotLayerHierarchy::updateAncestors(int const layerIndex, unsigned int const pivotIndex,
                                          tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap) {
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = ancestorIndicesMap.begin(); it1 != ancestorIndicesMap.end(); it1++) {
        int const ancestorLayerIndex = (*it1).first;
        tsl::sparse_set<unsigned int> const& ancestors = (*it1).second;

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = ancestors.begin(); it2 != ancestors.end(); it2++) {
            unsigned int const ancestorIndex = (*it2);

            // add ancestor and descendant
            (*_pivotLayers[layerIndex]).addAncestor(pivotIndex, ancestorLayerIndex, ancestorIndex);
            (*_pivotLayers[ancestorLayerIndex]).addDescendant(ancestorIndex, layerIndex, pivotIndex);

            // check to update dmax
            checkToUpdateDMax(layerIndex, pivotIndex, ancestorLayerIndex, ancestorIndex);
        }
    }
}

/**
 * @brief Add all descendants to the pivot. update dmax and umaxDescendants for newly added children.
 *
 * @param layerIndex
 * @param pivotIndex
 * @param descendantIndicesMap
 */
void PivotLayerHierarchy::updateDescendants(int const layerIndex, unsigned int const pivotIndex,
                                            tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& descendantIndicesMap) {
    // consider all layers below layerIndex
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = descendantIndicesMap.begin(); it1 != descendantIndicesMap.end(); it1++) {
        int const descendantLayerIndex = (*it1).first;
        tsl::sparse_set<unsigned int> const& descendants = (*it1).second;

        // go through each descendant and update data structures
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = descendants.begin(); it2 != descendants.end(); it2++) {
            unsigned int const descendantIndex = (*it2);

            // add ancestor and descendant
            (*_pivotLayers[layerIndex]).addDescendant(pivotIndex, descendantLayerIndex, descendantIndex);
            (*_pivotLayers[descendantLayerIndex]).addAncestor(descendantIndex, layerIndex, pivotIndex);

            // check to update dmax
            checkToUpdateDMax(descendantLayerIndex,descendantIndex,layerIndex,pivotIndex);
        }

        // should have added all descendants of pivotIndex in descendant layer index. 
        // now we have to compute umaxDescendants and vmaxDescendants of pivotIndex. 
        _computeUmaxDescendantMaxLinkDistance(layerIndex,pivotIndex,descendantLayerIndex);
        
        if (layerIndex == 0) {
            IndexPairDistanceStruct const& umaxDescendant = (*_pivotLayers[layerIndex]).get_umaxDescendantMaxLinkDistance(pivotIndex,descendantLayerIndex);
            IndexPairDistanceStruct const& umaxMaxLayer = (*_pivotLayers[layerIndex]).get_umaxMaximumOfDescendantMaxLinkDistance(descendantLayerIndex);

            if (umaxMaxLayer.isInvalid()) {
                (*_pivotLayers[0]).set_umaxMaximumOfDescendantMaxLinkDistance(descendantLayerIndex, umaxDescendant);
            } else if (umaxDescendant.distance > umaxMaxLayer.distance) {
                (*_pivotLayers[0]).set_umaxMaximumOfDescendantMaxLinkDistance(descendantLayerIndex, umaxDescendant);
            }
        }

        // current state umax,dmax (probably uninitialized since this is used for new queries)
        _computeVmaxDescendantMaxCoarseLinkDistance(layerIndex,pivotIndex,descendantLayerIndex);//(layerIndex,pivotIndex,descendantLayerIndex);

        if (layerIndex == 0) {
            IndexPairDistanceStruct const& vmaxDescendant = (*_pivotLayers[layerIndex]).get_vmaxDescendantMaxCoarseLinkDistance(pivotIndex,descendantLayerIndex);
            IndexPairDistanceStruct const& vmaxMaxLayer = (*_pivotLayers[layerIndex]).get_vmaxMaximumOfDescendantMaxCoarseLinkDistance(descendantLayerIndex);

            if (vmaxMaxLayer.isInvalid()) {
                (*_pivotLayers[0]).set_vmaxMaximumOfDescendantMaxCoarseLinkDistance(descendantLayerIndex, vmaxDescendant);
            } else if (vmaxDescendant.distance > vmaxMaxLayer.distance) {
                (*_pivotLayers[0]).set_vmaxMaximumOfDescendantMaxCoarseLinkDistance(descendantLayerIndex, vmaxDescendant);
            }
        }
    }

    return;
}

/**
 * @brief Check if new child pivotIndex of ancestorIndex requires dmax update
 *
 * @param layerIndex
 * @param pivotIndex
 * @param ancestorLayerIndex
 * @param ancestorIndex
 */
void PivotLayerHierarchy::checkToUpdateDMax(int const layerIndex, unsigned int const pivotIndex, int const ancestorLayerIndex,
                                            unsigned int const ancestorIndex) {
    // update dmaxMaxChildDistance_ of ancestor
    float const distance = getDistance(pivotIndex, ancestorIndex);
    if (distance > (*_pivotLayers[ancestorLayerIndex]).get_dmaxMaxChildDistance(ancestorIndex, layerIndex)) {
        (*_pivotLayers[ancestorLayerIndex]).set_dmaxMaxChildDistance(ancestorIndex, layerIndex, distance);

        // update dmaxMaximumOfMaxChildDistance_ of layer
        if (distance > (*_pivotLayers[ancestorLayerIndex]).get_dmaxMaximumOfMaxChildDistance(layerIndex)) {
            (*_pivotLayers[ancestorLayerIndex]).set_dmaxMaximumOfMaxChildDistance(layerIndex, distance);
        }
    }
}



/**
 * @brief Add given set of coarse neighbors of pivotIndex in fine-layer layerIndex
 *
 * @param layerIndex
 * @param pivotIndex
 * @param coarsePivotLayerNeighbors
 */
void PivotLayerHierarchy::addCoarsePivotLayerNeighbors(int const layerIndex, unsigned int const pivotIndex,
                                                       tsl::sparse_set<unsigned int> const& coarsePivotLayerNeighbors) {
    // add each coarse neighbor of pivotIndex
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarsePivotLayerNeighbors.begin(); it1 != coarsePivotLayerNeighbors.end(); it1++) {
        unsigned int coarsePivotIndex = (*it1);  // the coarse neighbor in layerIndex-1
        (*_pivotLayers[layerIndex]).add_coarsePivotLayerNeighbor(pivotIndex, coarsePivotIndex);

        // update pivotIndex's vmaxMaxCoarseLinkDistance
        _updateVmaxMaxCoarseLinkDistance(layerIndex, pivotIndex, coarsePivotIndex);
    }
}

/**
 * @brief add Q as coarse neighbors of pivots in finerPivotLayerNeighbors in layer layerIndex+1
 *
 * @param layerIndex
 * @param pivotIndex
 * @param finerPivotLayerNeighbors
 */
void PivotLayerHierarchy::addFinerCoarsePivotLayerNeighbors(int const layerIndex, unsigned int const pivotIndex,
                                                            tsl::sparse_set<unsigned int> const& finerPivotLayerNeighbors) {
    int finerLayerIndex = layerIndex + 1;

    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = finerPivotLayerNeighbors.begin(); it1 != finerPivotLayerNeighbors.end(); it1++) {
        unsigned int finerPivotIndex = (*it1);

        (*_pivotLayers[finerLayerIndex]).add_coarsePivotLayerNeighbor(finerPivotIndex, pivotIndex);

        // update finerPivotIndex's max coarse link distance, and propogate to ancestors
        _updateVmaxMaxCoarseLinkDistance(finerLayerIndex, finerPivotIndex, pivotIndex);
    }

    return;
}

/**
 * @brief Q in Layer layerIndex invalidated coarse neighbor links of fineIndex in layerIndex+1. remove.
 *
 * @param layerIndex
 * @param invalidatedPivotLinks
 */
void PivotLayerHierarchy::removeFinerCoarsePivotLayerNeighbors(
    int const layerIndex, tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& invalidatedPivotLinks) {
    int finerLayerIndex = layerIndex + 1;

    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = invalidatedPivotLinks.begin(); it1 != invalidatedPivotLinks.end(); it1++) {
        unsigned int const finePivotIndex = (*it1).first;
        tsl::sparse_set<unsigned int> const& invalidatedCoarseNeighbors = (*it1).second;
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = invalidatedCoarseNeighbors.begin(); it2 != invalidatedCoarseNeighbors.end(); it2++) {
            unsigned int invalidatedCoarseNeighborIndex = (*it2);
            (*_pivotLayers[finerLayerIndex]).remove_coarsePivotLayerNeighbor(finePivotIndex, invalidatedCoarseNeighborIndex);

            // this handles all vmax: maxCoarse,descendants,layer
            _checkToRecomputeVmaxMaxCoarseLinkDistance(finerLayerIndex, finePivotIndex, invalidatedCoarseNeighborIndex);
        }
    }

    return;
};



/**
 * @brief Get stats for each layer: ave degree, min, max, num links
 *
 * @param numLinks
 * @param aveDegree
 * @param minDegree
 * @param maxDegree
 */
void PivotLayerHierarchy::getDegreeStats(std::vector<unsigned int>& numLinks, std::vector<float>& aveDegree,
                                              std::vector<unsigned int>& minDegree, std::vector<unsigned int>& maxDegree) {
    numLinks.clear();
    aveDegree.clear();
    minDegree.clear();
    maxDegree.clear();
    numLinks.resize(_numberOfPivotLayers);
    aveDegree.resize(_numberOfPivotLayers);
    minDegree.resize(_numberOfPivotLayers);
    maxDegree.resize(_numberOfPivotLayers);

    for (int i = 0; i < _numberOfPivotLayers; i++) {
        unsigned int linkCount = 0;
        unsigned int min_degree = 100000000;
        unsigned int max_degree = 0;

        // go through each pivot in layer. find neighbors.
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> neighbors = *((*_pivotLayers[i]).get_pivotLayerNeighbors_ptr());
        unsigned int const numberOfPivots = (*_pivotLayers[i]).get_numberOfPivots();
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>::const_iterator it1;
        for (it1 = neighbors.begin(); it1 != neighbors.end(); it1++) {
            unsigned int lc = ((*it1).second).size();

            linkCount += lc;
            if (lc < min_degree) {
                min_degree = lc;
            }
            if (lc > max_degree) {
                max_degree = lc;
            }
        }

        numLinks[i] = ((unsigned int)linkCount / 2);  // counting each link twice
        aveDegree[i] = ((float)linkCount / numberOfPivots);
        minDegree[i] = min_degree;
        maxDegree[i] = max_degree;
    }

    return;
}