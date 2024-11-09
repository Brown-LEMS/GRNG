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
// 01-26-2022
#include "pivot-layer-hierarchy.hpp"

/**
 * @brief pivotIndex in layerIndex has new coarse neighbor, coarseIndex. update vmax of pivotIndex
 *
 * @param layerIndex
 * @param pivotIndex
 * @param coarsePivotIndex
 */
void PivotLayerHierarchy::_updateVmaxMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex, unsigned int const coarsePivotIndex) {
    float const pivotRadius = (*_pivotLayers[layerIndex]).get_pivotRadius(pivotIndex);
    float const coarseRadius = (*_pivotLayers[layerIndex - 1]).get_pivotRadius(coarsePivotIndex);
    float const distance12 = getDistance(pivotIndex, coarsePivotIndex);
    float const vmax_distance = _luneRadius(distance12, pivotRadius, coarseRadius);
    IndexPairDistanceStruct temp_ipds(pivotIndex, coarsePivotIndex, vmax_distance);

    // update vmaxMaxCoarseLinkDistance_ of pivotIndex1Radius
    if ((*_pivotLayers[layerIndex]).get_vmaxMaxCoarseLinkDistance(pivotIndex).isInvalid()) {
        (*_pivotLayers[layerIndex]).set_vmaxMaxCoarseLinkDistance(pivotIndex, temp_ipds);
    } else if (vmax_distance > (*_pivotLayers[layerIndex]).get_vmaxMaxCoarseLinkDistance(pivotIndex).distance) {
        (*_pivotLayers[layerIndex]).set_vmaxMaxCoarseLinkDistance(pivotIndex, temp_ipds);
    }

    _updateVmaxDescendantMaxCoarseLinkDistance(layerIndex, pivotIndex);
}

/**
 * @brief pivotIndex in layerIndex has new coarse neighbor. update vmax for all ancestors of x_i
 *
 * @param layerIndex
 * @param pivotIndex
 */
void PivotLayerHierarchy::_updateVmaxDescendantMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex) {
    //if (layerIndex == 0) return;  // no pivots above this layer

    // get vmax of pivotIndex, max distance to a coarse neighbor
    IndexPairDistanceStruct const& vmax_pivotIndex = (*_pivotLayers[layerIndex]).get_vmaxMaxCoarseLinkDistance(pivotIndex);

    // consider each ancestor of pivotIndex
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap =
        (*_pivotLayers[layerIndex]).get_ancestorPivotIndicesMap(pivotIndex);
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = ancestorIndicesMap.begin(); it1 != ancestorIndicesMap.end(); it1++) {
        int const ancestorLayerIndex = (*it1).first;
        tsl::sparse_set<unsigned int> const& ancestors = (*it1).second;

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = ancestors.begin(); it2 != ancestors.end(); it2++) {
            unsigned int const ancestorIndex = (*it2);

            // current value of vmax descendant
            IndexPairDistanceStruct const& vmax_ancestor_current =
                (*_pivotLayers[ancestorLayerIndex]).get_vmaxDescendantMaxCoarseLinkDistance(ancestorIndex, layerIndex);

            // compute new value of vmax descendant
            float const vmax_descendant_value = getDistance(pivotIndex, ancestorIndex) + vmax_pivotIndex.distance;
            IndexPairDistanceStruct potential_vmax_descendant(vmax_pivotIndex.index1, vmax_pivotIndex.index2, vmax_descendant_value);

            // check if greater
            bool flag_update_max = false;
            if (vmax_ancestor_current.isInvalid()) {
                (*_pivotLayers[ancestorLayerIndex]).set_vmaxDescendantMaxCoarseLinkDistance(ancestorIndex, layerIndex, potential_vmax_descendant);
                flag_update_max = true;
            } else if (vmax_descendant_value > vmax_ancestor_current.distance) {
                (*_pivotLayers[ancestorLayerIndex]).set_vmaxDescendantMaxCoarseLinkDistance(ancestorIndex, layerIndex, potential_vmax_descendant);
                flag_update_max = true;
            }

            // update vmaxDescendantMaxCoarseLinkDistance if necessary
            if (ancestorLayerIndex == 0 && flag_update_max) {
                IndexPairDistanceStruct const& vmax_max_current = (*_pivotLayers[0]).get_vmaxMaximumOfDescendantMaxCoarseLinkDistance(layerIndex);

                if (vmax_max_current.isInvalid()) {
                    (*_pivotLayers[0]).set_vmaxMaximumOfDescendantMaxCoarseLinkDistance(layerIndex, potential_vmax_descendant);
                } else if (vmax_descendant_value > vmax_max_current.distance) {
                    (*_pivotLayers[0]).set_vmaxMaximumOfDescendantMaxCoarseLinkDistance(layerIndex, potential_vmax_descendant);
                }
            }
        }
    }
}

/**
 * @brief coarse neighbor of pivotIndex in layerIndex was removed. check if recompute vmax for pivotIndex.
 *
 * @param layerIndex
 * @param pivotIndex
 * @param coarseIndex
 */
void PivotLayerHierarchy::_checkToRecomputeVmaxMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex,
                                                                     unsigned int const coarseIndex) {
    if ((*_pivotLayers[layerIndex]).get_vmaxMaxCoarseLinkDistance(pivotIndex).isExactMatch(pivotIndex, coarseIndex)) {
        _computeVmaxMaxCoarseLinkDistance(layerIndex, pivotIndex);
    }

    // coarse link pivotIndex -> coarseIndex was eliminated. does this impact ancestors?
    checkAncestorsToRecomputeVmaxDescendantCoarse(layerIndex, pivotIndex, coarseIndex);

    // check this too!!
    _checkToRecomputeVmaxMaximumOfDescendantCoarse(layerIndex,pivotIndex,coarseIndex);
}



/**
 * @brief finePivot in layerIndex no longer has coarse neighbor coarsePivot. check ancestors of finePivot to recompute.
 *
 * @param layerIndex
 * @param finePivot
 * @param coarsePivot
 */
void PivotLayerHierarchy::checkAncestorsToRecomputeVmaxDescendantCoarse(int const layerIndex, unsigned int const finePivot,
                                                                        unsigned int const coarsePivot) {
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap =
        (*_pivotLayers[layerIndex]).get_ancestorPivotIndicesMap(finePivot);

    tsl::sparse_map<int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = ancestorIndicesMap.begin(); it1 != ancestorIndicesMap.end(); it1++) {
        unsigned int const coarserLayerIndex = (*it1).first;
        tsl::sparse_set<unsigned int> const& ancestorIndices = (*it1).second;

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = ancestorIndices.begin(); it2 != ancestorIndices.end(); it2++) {
            unsigned int const ancestorIndex = (*it2);
            IndexPairDistanceStruct const& ancestorVmaxDescendant =
                (*_pivotLayers[coarserLayerIndex]).get_vmaxDescendantMaxCoarseLinkDistance(ancestorIndex, layerIndex);
            if (ancestorVmaxDescendant.isExactMatch(finePivot, coarsePivot)) {
                _computeVmaxDescendantMaxCoarseLinkDistance(coarserLayerIndex, ancestorIndex, layerIndex);
            }
        }
    }
}

void PivotLayerHierarchy::_checkToRecomputeVmaxMaximumOfDescendantCoarse(int const layerIndex, unsigned int const fineIndex, unsigned int const coarseIndex) {
    IndexPairDistanceStruct const vmaxMaxLayer = (*_pivotLayers[0]).get_vmaxMaximumOfDescendantMaxCoarseLinkDistance(layerIndex);

    if (vmaxMaxLayer.isExactMatch(fineIndex, coarseIndex)) {
        IndexPairDistanceStruct vmax_maxValue;  // max val holder

        // search through all existing vmaxDescendantMaxCoarseLinkDistance_ for max
        tsl::sparse_set<unsigned int> const& topLayerPivotIndices = *(*_pivotLayers[0]).get_pivotIndices_ptr();
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = topLayerPivotIndices.begin(); it1 != topLayerPivotIndices.end(); it1++) {
            unsigned int const topLayerPivotIndex = (*it1);
            IndexPairDistanceStruct const& topLayerPivotVmaxValue =
                (*_pivotLayers[0]).get_vmaxDescendantMaxCoarseLinkDistance(topLayerPivotIndex, layerIndex);

            if (!topLayerPivotVmaxValue.isInvalid()) {
                if (topLayerPivotVmaxValue.distance > vmax_maxValue.distance) {
                    vmax_maxValue = topLayerPivotVmaxValue;
                }
            }
        }
        (*_pivotLayers[0]).set_vmaxMaximumOfDescendantMaxCoarseLinkDistance(layerIndex, vmax_maxValue);
    }
}


/**
 * @brief Compute max coarse link for pivotIndex in layerIndex
 *
 * @param layerIndex
 * @param pivotIndex
 */
void PivotLayerHierarchy::_computeVmaxMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex) {
    float const pivotRadius = (*_pivotLayers[layerIndex]).get_pivotRadius(pivotIndex);

    // reinitialize vmaxMaxCoarseLinkDistance_
    IndexPairDistanceStruct vmax_maxValue;

    // Find max link distance of all neighbors of index1
    tsl::sparse_set<unsigned int> const& coarseNeighborIndices = (*_pivotLayers[layerIndex]).get_coarsePivotLayerNeighbors(pivotIndex);
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = coarseNeighborIndices.begin(); it1 != coarseNeighborIndices.end(); it1++) {
        unsigned int const coarseIndex = (*it1);
        float const coarseRadius = (*_pivotLayers[layerIndex - 1]).get_pivotRadius(coarseIndex);

        // compute the link distance for this neighbor
        float const distance12 = getDistance(pivotIndex, coarseIndex);
        float const vmax_distance = _luneRadius(distance12, pivotRadius, coarseRadius);

        // update vmax_maxValue if the value is greater
        if (vmax_maxValue.isInvalid()) {
            vmax_maxValue = IndexPairDistanceStruct(pivotIndex, coarseIndex, vmax_distance);
        } else if (vmax_distance > vmax_maxValue.distance) {
            vmax_maxValue = IndexPairDistanceStruct(pivotIndex, coarseIndex, vmax_distance);
        }
    }

    // update vmaxMaxCoarseLinkDistance_ with max value for neighbor
    (*_pivotLayers[layerIndex]).set_vmaxMaxCoarseLinkDistance(pivotIndex, vmax_maxValue);
}


/**
 * @brief Recompute vmaxDescendant of pivotIndex in layer index. refer to descendants in finerLayerIndex
 *
 * @param layerIndex
 * @param ancestorLayerIndex
 * @param ancestorIndex
 */
void PivotLayerHierarchy::_computeVmaxDescendantMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex,
                                                                      int const finerLayerIndex) {
    IndexPairDistanceStruct vmaxDescendant_maxValue;

    tsl::sparse_set<unsigned int> const& descendants = (*_pivotLayers[layerIndex]).get_descendantPivotIndices(pivotIndex, finerLayerIndex);

    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = descendants.begin(); it1 != descendants.end(); it1++) {
        unsigned int const descendantIndex = (*it1);

        // compute the vmaxDescendantMaxCoarseLinkDistance for this relationship
        IndexPairDistanceStruct const& vmax_descendantIndex = (*_pivotLayers[finerLayerIndex]).get_vmaxMaxCoarseLinkDistance(descendantIndex);
        float const vmax_descendant_distance = getDistance(pivotIndex, descendantIndex) + vmax_descendantIndex.distance;
        IndexPairDistanceStruct vmax_descendant_temp(vmax_descendantIndex.index1, vmax_descendantIndex.index2, vmax_descendant_distance);

        // update vmax_maxValue if applicable
        if (vmaxDescendant_maxValue.isInvalid()) {
            vmaxDescendant_maxValue = vmax_descendant_temp;
        } else if (vmax_descendant_distance > vmaxDescendant_maxValue.distance) {
            vmaxDescendant_maxValue = vmax_descendant_temp;
        }
    }

    (*_pivotLayers[layerIndex]).set_vmaxDescendantMaxCoarseLinkDistance(pivotIndex, finerLayerIndex, vmaxDescendant_maxValue);
}



