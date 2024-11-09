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
#include "pivot-layer.hpp"

PivotLayer::PivotLayer(int const pivotLayerID, int const numberOfPivotLayers,
                       std::shared_ptr<SparseMatrix> const SparseMatrix,
                       std::string luneType, bool verbose) : _pivotLayerID(pivotLayerID), _numberOfPivotLayers(numberOfPivotLayers), _sparseMatrix(SparseMatrix), _luneType(luneType), _verbose(verbose) {
    // initialize layer data structures
    (*_dmaxMaximumOfMaxChildDistance).clear();
    for (int finerLayerIndex = _pivotLayerID + 1; finerLayerIndex < numberOfPivotLayers; finerLayerIndex++) {
        (*_dmaxMaximumOfMaxChildDistance).emplace(finerLayerIndex, 0.0f);
    }
    if (_pivotLayerID == 0) {
        (*_umaxMaximumOfDescendantMaxLinkDistance).clear();
        for (int finerLayerIndex = _pivotLayerID + 1; finerLayerIndex < numberOfPivotLayers; finerLayerIndex++) {
            (*_umaxMaximumOfDescendantMaxLinkDistance).emplace(finerLayerIndex, IndexPairDistanceStruct());
        }
        (*_vmaxMaximumOfDescendantMaxCoarseLinkDistance).clear();
        for (int finerLayerIndex = _pivotLayerID + 1; finerLayerIndex < numberOfPivotLayers; finerLayerIndex++) {
            (*_vmaxMaximumOfDescendantMaxCoarseLinkDistance).emplace(finerLayerIndex, IndexPairDistanceStruct());
        }
    }

    // set lune type (defaults to Berk)
    if (_luneType == "Cole") {
        _luneIndex = 1;
    }
};

// copy constructor, copy all data structures
PivotLayer::PivotLayer(PivotLayer const& layer) : _pivotLayerID(layer.get_pivotLayerID()), _numberOfPivotLayers(layer.get_numberOfPivotLayers()), _sparseMatrix(layer.get_sparseMatrix()), _luneType(layer.get_luneType()) {
    // set lune type
    if (_luneType == "Cole") {
        _luneIndex = 1;
    };

    // copy all data structures
    _pivotIndices = layer.get_pivotIndices_ptr();
    _pivotRadii = layer.get_pivotRadii_ptr();
    _pivotLayerNeighbors = layer.get_pivotLayerNeighbors_ptr();
    _coarsePivotLayerNeighbors = layer.get_coarsePivotLayerNeighbors_ptr();
    _ancestorPivotIndices = layer.get_ancestorPivotIndices_ptr();
    _descendantPivotIndices = layer.get_descendantPivotIndices_ptr();
    _dmaxMaxChildDistance = layer.get_dmaxMaxChildDistance_ptr();
    _dmaxMaximumOfMaxChildDistance = layer.get_dmaxMaximumOfMaxChildDistance_ptr();
    _umaxMaxLinkDistance = layer.get_umaxMaxLinkDistance_ptr();
    _umaxDescendantMaxLinkDistance = layer.get_umaxDescendantMaxLinkDistance_ptr();
    _umaxMaximumOfDescendantMaxLinkDistance = layer.get_umaxMaximumOfDescendantMaxLinkDistance_ptr();
    _vmaxMaxCoarseLinkDistance = layer.get_vmaxMaxCoarseLinkDistance_ptr();
    _vmaxDescendantMaxCoarseLinkDistance = layer.get_vmaxDescendantMaxCoarseLinkDistance_ptr();
    _vmaxMaximumOfDescendantMaxCoarseLinkDistance = layer.get_vmaxMaximumOfDescendantMaxCoarseLinkDistance_ptr();
}

void PivotLayer::addPivotToLayer(unsigned int const pivotIndex, float const radius) {
    (*_pivotIndices).insert(pivotIndex);
    (*_pivotRadii).emplace(pivotIndex, radius);
}

void PivotLayer::initializePivotInLayer(unsigned int const pivotIndex) {
    (*_pivotLayerNeighbors).emplace(pivotIndex, tsl::sparse_set<unsigned int>{});

    // top layer pivots
    if (_pivotLayerID == 0) {
        // lower layer structures
        (*_descendantPivotIndices).emplace(pivotIndex, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{});
        (*_dmaxMaxChildDistance).emplace(pivotIndex, tsl::sparse_map<int, float>{});
        (*_umaxDescendantMaxLinkDistance).emplace(pivotIndex, tsl::sparse_map<int, IndexPairDistanceStruct>{});
        (*_vmaxDescendantMaxCoarseLinkDistance).emplace(pivotIndex, tsl::sparse_map<int, IndexPairDistanceStruct>{});
        for (int finerLayerIndex = _pivotLayerID + 1; finerLayerIndex < _numberOfPivotLayers; finerLayerIndex++) {
            (*_descendantPivotIndices)[pivotIndex].emplace(finerLayerIndex, tsl::sparse_set<unsigned int>{});
            (*_dmaxMaxChildDistance)[pivotIndex].emplace(finerLayerIndex, 0.0f);
            (*_umaxDescendantMaxLinkDistance)[pivotIndex].emplace(finerLayerIndex, IndexPairDistanceStruct());
            (*_vmaxDescendantMaxCoarseLinkDistance)[pivotIndex].emplace(finerLayerIndex, IndexPairDistanceStruct());
        }
    } else if (_pivotLayerID == _numberOfPivotLayers - 1) {
        // current layer structures
        (*_coarsePivotLayerNeighbors).emplace(pivotIndex, tsl::sparse_set<unsigned int>{});
        (*_umaxMaxLinkDistance).emplace(pivotIndex, IndexPairDistanceStruct());
        (*_vmaxMaxCoarseLinkDistance).emplace(pivotIndex, IndexPairDistanceStruct());

        // upper layer structures
        (*_ancestorPivotIndices).emplace(pivotIndex, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{});
        for (int coarserLayerIndex = 0; coarserLayerIndex < _pivotLayerID; coarserLayerIndex++) {
            (*_ancestorPivotIndices)[pivotIndex].emplace(coarserLayerIndex, tsl::sparse_set<unsigned int>{});
        }
    } else {  // middle layer pivots

        // current layer structures
        (*_coarsePivotLayerNeighbors).emplace(pivotIndex, tsl::sparse_set<unsigned int>{});
        (*_umaxMaxLinkDistance).emplace(pivotIndex, IndexPairDistanceStruct());
        (*_vmaxMaxCoarseLinkDistance).emplace(pivotIndex, IndexPairDistanceStruct());

        // upper layer structures
        (*_ancestorPivotIndices).emplace(pivotIndex, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{});
        for (int coarserLayerIndex = 0; coarserLayerIndex < _pivotLayerID; coarserLayerIndex++) {
            (*_ancestorPivotIndices)[pivotIndex].emplace(coarserLayerIndex, tsl::sparse_set<unsigned int>{});
        }

        // lower layer structures
        (*_descendantPivotIndices).emplace(pivotIndex, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>{});
        (*_dmaxMaxChildDistance).emplace(pivotIndex, tsl::sparse_map<int, float>{});
        (*_umaxDescendantMaxLinkDistance).emplace(pivotIndex, tsl::sparse_map<int, IndexPairDistanceStruct>{});
        (*_vmaxDescendantMaxCoarseLinkDistance).emplace(pivotIndex, tsl::sparse_map<int, IndexPairDistanceStruct>{});
        for (int finerLayerIndex = _pivotLayerID + 1; finerLayerIndex < _numberOfPivotLayers; finerLayerIndex++) {
            (*_descendantPivotIndices)[pivotIndex].emplace(finerLayerIndex, tsl::sparse_set<unsigned int>{});
            (*_dmaxMaxChildDistance)[pivotIndex].emplace(finerLayerIndex, 0.0f);
            (*_umaxDescendantMaxLinkDistance)[pivotIndex].emplace(finerLayerIndex, IndexPairDistanceStruct());
            (*_vmaxDescendantMaxCoarseLinkDistance)[pivotIndex].emplace(finerLayerIndex, IndexPairDistanceStruct());
        }
    }
}

// add edge to pivotLayerNeighbors and update umax
void PivotLayer::addLinkToLayer(unsigned int const pivotIndex1, unsigned int const pivotIndex2) {
    (*_pivotLayerNeighbors)[pivotIndex1].insert(pivotIndex2);
    (*_pivotLayerNeighbors)[pivotIndex2].insert(pivotIndex1);

    if (_pivotLayerID > 0) {
        _updateUmaxMaxLinkDistance(pivotIndex1, pivotIndex2);
    }
}

// remove edge from pivotLayerNeighbors and recompute umax if necessary
void PivotLayer::removeLinkFromLayer(unsigned int const pivotIndex1, unsigned int const pivotIndex2) {
    (*_pivotLayerNeighbors)[pivotIndex1].erase(pivotIndex2);
    (*_pivotLayerNeighbors)[pivotIndex2].erase(pivotIndex1);

    if (_pivotLayerID > 0) {
        if ((*_umaxMaxLinkDistance)[pivotIndex1].isLinkMatch(pivotIndex1, pivotIndex2)) {
            _computeUmaxMaxLinkDistance(pivotIndex1);
        }
        if ((*_umaxMaxLinkDistance)[pivotIndex2].isLinkMatch(pivotIndex2, pivotIndex1)) {
            _computeUmaxMaxLinkDistance(pivotIndex2);
        }
    }
}

/**
 * @brief add a new coarse neighbor, update vMax 
 * 
 * @param pivotIndex 
 * @param coarseNeighbor 
 */
void PivotLayer::add_coarsePivotLayerNeighbor(unsigned int const pivotIndex, unsigned int const coarseNeighbor) {
    (*_coarsePivotLayerNeighbors)[pivotIndex].insert(coarseNeighbor);
}

/**
 * @brief remove a coarse neighbor, recompute vMax if necessary
 * 
 * @param pivotIndex 
 * @param coarseNeighbor 
 */
void PivotLayer::remove_coarsePivotLayerNeighbor(unsigned int const pivotIndex, unsigned int const coarseNeighbor) {
    (*_coarsePivotLayerNeighbors)[pivotIndex].erase(coarseNeighbor);
}


// update umax for each index of the new link
void PivotLayer::_updateUmaxMaxLinkDistance(unsigned int const pivotIndex1, unsigned int const pivotIndex2) {
    float const pivotIndex1Radius = (*_pivotRadii)[pivotIndex1];
    float const pivotIndex2Radius = (*_pivotRadii)[pivotIndex2];
    float const distance12 = getDistance(pivotIndex1, pivotIndex2);

    // update umaxMaxLinkDistance_ of pivotIndex1Radius
    float const umax_distance1 = _luneRadius(distance12, pivotIndex1Radius, pivotIndex2Radius);
    if ((*_umaxMaxLinkDistance)[pivotIndex1].isInvalid()) {
        (*_umaxMaxLinkDistance)[pivotIndex1] = IndexPairDistanceStruct(pivotIndex1, pivotIndex2, umax_distance1);
    } else if (umax_distance1 > (*_umaxMaxLinkDistance)[pivotIndex1].distance) {
        (*_umaxMaxLinkDistance)[pivotIndex1] = IndexPairDistanceStruct(pivotIndex1, pivotIndex2, umax_distance1);
    }

    // update umaxMaxLinkDistance_ of index2
    float const umax_distance2 = _luneRadius(distance12, pivotIndex2Radius, pivotIndex1Radius);
    if ((*_umaxMaxLinkDistance)[pivotIndex2].isInvalid()) {
        (*_umaxMaxLinkDistance)[pivotIndex2] = IndexPairDistanceStruct(pivotIndex2, pivotIndex1, umax_distance2);
    } else if (umax_distance2 > (*_umaxMaxLinkDistance)[pivotIndex2].distance) {
        (*_umaxMaxLinkDistance)[pivotIndex2] = IndexPairDistanceStruct(pivotIndex2, pivotIndex1, umax_distance2);
    }
}

// compute umax for index by looking at all current links
void PivotLayer::_computeUmaxMaxLinkDistance(unsigned int const pivotIndex1) {
    float const pivotIndex1Radius = (*_pivotRadii)[pivotIndex1];

    // reinitialize umaxMaxLinkDistance_
    (*_umaxMaxLinkDistance)[pivotIndex1] = IndexPairDistanceStruct();
    IndexPairDistanceStruct umax_maxValue;

    // Find max link distance of all neighbors of index1
    tsl::sparse_set<unsigned int> const& pivotIndex1NeighborIndices = (*_pivotLayerNeighbors)[pivotIndex1];
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndex1NeighborIndices.begin(); it1 != pivotIndex1NeighborIndices.end(); it1++) {
        unsigned int const pivotIndex2 = (*it1);
        float const pivotIndex2Radius = (*_pivotRadii)[pivotIndex2];

        // compute the link distance for this neighbor
        float const distance12 = getDistance(pivotIndex1, pivotIndex2);
        float const umax_distance = _luneRadius(distance12, pivotIndex1Radius, pivotIndex2Radius);

        // update umax_maxValue if the value is greater
        if (umax_maxValue.isInvalid()) {
            umax_maxValue = IndexPairDistanceStruct(pivotIndex1, pivotIndex2, umax_distance);
        } else if (umax_distance > umax_maxValue.distance) {
            umax_maxValue = IndexPairDistanceStruct(pivotIndex1, pivotIndex2, umax_distance);
        }
    }

    // update umaxMaxLinkDistance_ with max value for neighbor
    (*_umaxMaxLinkDistance)[pivotIndex1] = umax_maxValue;
}

void PivotLayer::computeBruteForcePivotLayerGraph() {
    // sort indices as to not repeat checks for neighbors (if index2 < index1)
    std::vector<unsigned int> pivotIndices((*_pivotIndices).begin(), (*_pivotIndices).end());
    std::sort(pivotIndices.begin(), pivotIndices.end());

    // consider each index1
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex1 = (*it1);
        float const pivotIndex1Radius = (*_pivotRadii)[pivotIndex1];

        // consider each index2 as a neighbor of pivotIndex1
        std::vector<unsigned int>::const_iterator it2;
        for (it2 = pivotIndices.begin(); it2 != pivotIndices.end(); it2++) {
            unsigned int const pivotIndex2 = (*it2);
            float const pivotIndex2Radius = (*_pivotRadii)[pivotIndex2];

            if (pivotIndex1 >= pivotIndex2) continue;  // this is an issue..... fuckk???

            bool flag_isNeighbor = true;
            float const distance12 = getDistance(pivotIndex1, pivotIndex2);

            // if no lune, no interference --> automatic neighbor
            if (_luneRadius(distance12, pivotIndex1Radius, pivotIndex2Radius) <= 0 ||
                _luneRadius(distance12, pivotIndex2Radius, pivotIndex1Radius) <= 0) {
                addLinkToLayer(pivotIndex1, pivotIndex2);
                continue;
            }

            // check for interference by other pivots
            std::vector<unsigned int>::const_iterator it3;
            for (it3 = pivotIndices.begin(); it3 != pivotIndices.end(); it3++) {
                unsigned int const pivotIndex3 = (*it3);
                if (pivotIndex3 == pivotIndex1 || pivotIndex3 == pivotIndex2) continue;

                // check both sides for interference
                float const distance13 = getDistance(pivotIndex1, pivotIndex3);
                if (distance13 < _luneRadius(distance12, pivotIndex1Radius, pivotIndex2Radius)) {
                    float const distance23 = getDistance(pivotIndex2, pivotIndex3);
                    if (distance23 < _luneRadius(distance12, pivotIndex2Radius, pivotIndex1Radius)) {
                        flag_isNeighbor = false;
                        break;
                    }
                }
            }

            if (flag_isNeighbor) {
                addLinkToLayer(pivotIndex1, pivotIndex2);
            }
        }
    }
}


