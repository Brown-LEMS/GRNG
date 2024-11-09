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
 * @brief update umaxDescendant based on newly added pivot pivotIndex
 *
 * @param layerIndex
 * @param pivotIndex
 */
void PivotLayerHierarchy::updateUmaxDescendantMaxLinkDistance_(int const layerIndex, unsigned int const pivotIndex) {
    IndexPairDistanceStruct const& umax_pivotIndex = (*_pivotLayers[layerIndex]).get_umaxMaxLinkDistance(pivotIndex);

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

            // current value of umax descendant
            IndexPairDistanceStruct const& umax_ancestor_current =
                (*_pivotLayers[ancestorLayerIndex]).get_umaxDescendantMaxLinkDistance(ancestorIndex, layerIndex);

            // compute new value of umax descendant
            float const umax_descendant_value = getDistance(pivotIndex, ancestorIndex) + umax_pivotIndex.distance;
            IndexPairDistanceStruct potential_umax_descendant(umax_pivotIndex.index1, umax_pivotIndex.index2, umax_descendant_value);

            // check if greater
            bool flag_update_max = false;
            if (umax_ancestor_current.isInvalid()) {
                (*_pivotLayers[ancestorLayerIndex]).set_umaxDescendantMaxLinkDistance(ancestorIndex, layerIndex, potential_umax_descendant);
                flag_update_max = true;
            } else if (umax_descendant_value > umax_ancestor_current.distance) {
                (*_pivotLayers[ancestorLayerIndex]).set_umaxDescendantMaxLinkDistance(ancestorIndex, layerIndex, potential_umax_descendant);
                flag_update_max = true;
            }

            // update umaxDescendantMaxLinkDistance if necessary
            if (ancestorLayerIndex == 0 && flag_update_max) {
                IndexPairDistanceStruct const& umax_max_current = (*_pivotLayers[0]).get_umaxMaximumOfDescendantMaxLinkDistance(layerIndex);

                if (umax_max_current.isInvalid()) {
                    (*_pivotLayers[0]).set_umaxMaximumOfDescendantMaxLinkDistance(layerIndex, potential_umax_descendant);
                } else if (umax_descendant_value > umax_max_current.distance) {
                    (*_pivotLayers[0]).set_umaxMaximumOfDescendantMaxLinkDistance(layerIndex, potential_umax_descendant);
                }
            }
        }
    }
}


/**
 * @brief Checking to update umaxDescendant of all ancestors after removed link
 *
 * @param layerIndex
 * @param pivotIndex1
 * @param pivotIndex2
 */
void PivotLayerHierarchy::checkAncestorsToRecomputeUmaxDescendant(int const layerIndex, unsigned int const pivotIndex1,
                                                                  unsigned int const pivotIndex2) {
    // go through each ancestor of pivotIndex1
    for (int coarserLayerIndex = 0; coarserLayerIndex < layerIndex; coarserLayerIndex++) {
        tsl::sparse_set<unsigned int> const& ancestorIndices = (*_pivotLayers[layerIndex]).get_ancestorPivotIndices(pivotIndex1, coarserLayerIndex);
        tsl::sparse_set<unsigned int>::const_iterator it3;
        for (it3 = ancestorIndices.begin(); it3 != ancestorIndices.end(); it3++) {
            unsigned int const ancestorIndex = (*it3);

            // check if umaxDescendantMaxLinkDistance impacted by removal of link (pivotIndex1,pivotIndex2)
            // if so, recompute umaxDescendantMaxLinkDistance entirely
            IndexPairDistanceStruct const& ancestorUmaxDescendant =
                (*_pivotLayers[coarserLayerIndex]).get_umaxDescendantMaxLinkDistance(ancestorIndex, layerIndex);
            if (ancestorUmaxDescendant.isLinkMatch(pivotIndex1, pivotIndex2)) {
                _computeUmaxDescendantMaxLinkDistance(coarserLayerIndex, ancestorIndex, layerIndex);
            }
        }
    }
}


/**
 * @brief recompute umaxDescendantMaxLinkDistance based on umaxLinkDistance recomputed via removed link
 *
 * @param layerIndex
 * @param ancestorLayerIndex
 * @param ancestorIndex
 */
void PivotLayerHierarchy::_computeUmaxDescendantMaxLinkDistance(int const layerIndex, unsigned int const ancestorIndex, int const finerLayerIndex) {
    // holder for max value
    IndexPairDistanceStruct umaxDescendant_maxValue;

    tsl::sparse_set<unsigned int> const& descendants = (*_pivotLayers[layerIndex]).get_descendantPivotIndices(ancestorIndex, finerLayerIndex);

    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = descendants.begin(); it1 != descendants.end(); it1++) {
        unsigned int const descendantIndex = (*it1);

        // compute the umaxDescendantMaxLinkDistance for this relationship
        IndexPairDistanceStruct const& umax_descendantIndex = (*_pivotLayers[finerLayerIndex]).get_umaxMaxLinkDistance(descendantIndex);
        float const umax_descendant_distance = getDistance(ancestorIndex, descendantIndex) + umax_descendantIndex.distance;
        IndexPairDistanceStruct umax_descendant_temp(umax_descendantIndex.index1, umax_descendantIndex.index2, umax_descendant_distance);

        // update umax_maxValue if applicable
        if (umaxDescendant_maxValue.isInvalid()) {
            umaxDescendant_maxValue = umax_descendant_temp;
        } else if (umax_descendant_distance > umaxDescendant_maxValue.distance) {
            umaxDescendant_maxValue = umax_descendant_temp;
        }
    }

    (*_pivotLayers[layerIndex]).set_umaxDescendantMaxLinkDistance(ancestorIndex, finerLayerIndex, umaxDescendant_maxValue);

    // check if need to update / compute umaxDescendant Max Link Distance
}



/**
 * @brief check if removed link impacts umaxDescendant of ancestors
 *
 * @param layerIndex
 * @param pivotIndex1
 * @param pivotIndex2
 */
void PivotLayerHierarchy::checkToRecomputeUmaxMaximumOfDescendant(int const layerIndex, unsigned int const pivotIndex1,
                                                                  unsigned int const pivotIndex2) {
    IndexPairDistanceStruct const umaxMaxLayer = (*_pivotLayers[0]).get_umaxMaximumOfDescendantMaxLinkDistance(layerIndex);
    if (umaxMaxLayer.isLinkMatch(pivotIndex1, pivotIndex2)) {
        IndexPairDistanceStruct umax_maxValue;  // max val holder

        // search through all existing umaxDescendantMaxLinkDistance_ for max
        tsl::sparse_set<unsigned int> const& topLayerPivotIndices = *(*_pivotLayers[0]).get_pivotIndices_ptr();
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = topLayerPivotIndices.begin(); it1 != topLayerPivotIndices.end(); it1++) {
            unsigned int const topLayerPivotIndex = (*it1);
            IndexPairDistanceStruct const& topLayerPivotUmaxValue =
                (*_pivotLayers[0]).get_umaxDescendantMaxLinkDistance(topLayerPivotIndex, layerIndex);

            if (!topLayerPivotUmaxValue.isInvalid()) {
                if (topLayerPivotUmaxValue.distance > umax_maxValue.distance) {
                    umax_maxValue = topLayerPivotUmaxValue;
                }
            }
        }
        (*_pivotLayers[0]).set_umaxMaximumOfDescendantMaxLinkDistance(layerIndex, umax_maxValue);
    }
}