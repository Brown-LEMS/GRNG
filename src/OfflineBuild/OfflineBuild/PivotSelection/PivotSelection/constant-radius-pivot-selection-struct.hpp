#ifndef ConstantRadiusPivotSelectionStruct_hpp
#define ConstantRadiusPivotSelectionStruct_hpp
/*
Copyright 2019, Brown University, Providence, RI.
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
// February 1st, 2022

#include <memory.h>
#include <tsl/sparse_map.h>
#include <tsl/sparse_set.h>

#include <vector>

/**
 * @brief Struct used during pivot selection to hold pivot indices in each layer and pivot domain members
 *
 */
class ConstantRadiusPivotSelectionStruct {
  public:
    ConstantRadiusPivotSelectionStruct(int const pivotLayerID, int const numberOfPivotLayers, float const radius) : pivotLayerID(pivotLayerID), numberOfPivotLayers(numberOfPivotLayers), pivotLayerRadius(radius){};
    ~ConstantRadiusPivotSelectionStruct(){};

    /**
     * @brief add pivot to the layer
     *
     * @param pivotIndex
     */
    inline void addPivot(unsigned int pivotIndex) {
        (*pivotIndices).insert(pivotIndex);
        if (pivotLayerID < numberOfPivotLayers - 2) {
            (*pivotMemberIndices)[pivotIndex] = tsl::sparse_set<unsigned int>{};
        }
    }

    /**
     * @brief add pivot as a member of parentIndex pivot domain
     *
     * @param parentIndex
     * @param pivotIndex
     */
    inline void addPivotMember(unsigned int parentIndex, unsigned int pivotIndex) {
        (*pivotMemberIndices)[parentIndex].insert(pivotIndex);
    }

    int const pivotLayerID = 0;
    int const numberOfPivotLayers = 1;
    float const pivotLayerRadius = 0.0f;

    /**
     * @brief set of all members of the pivot layer
     *
     */
    std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices =
        std::make_shared<tsl::sparse_set<unsigned int>>();

    /**
     * @brief map to contain pivot domain of each pivot in the layer
     *
     */
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> pivotMemberIndices =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>>();
};

#endif  // ConstantRadiusPivotSelectionStruct_hpp