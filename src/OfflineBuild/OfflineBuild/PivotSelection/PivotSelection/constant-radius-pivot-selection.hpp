#ifndef ConstantRadiusPivotSelection_hpp
#define ConstantRadiusPivotSelection_hpp
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

#include <vector>

#include "pivot-selection-struct.hpp"
#include "sparse-matrix.hpp"

/**
 * @brief Class used to perform Greedy Pivot Selection for the Hierarchical RNG
 *
 */
class ConstantRadiusPivotSelection {
  public:
    ConstantRadiusPivotSelection(unsigned int const datasetSize, std::vector<float> const& pivotRadiusVector,
                                 std::shared_ptr<SparseMatrix> const& sparseMatrix);
    ~ConstantRadiusPivotSelection(){};

    /**
     * @brief perform the greedy pivot selection in top-down fashion. return pivots belonging to each layer
     *
     * @return std::vector<PivotSelectionStruct>
     */
    std::vector<PivotSelectionStruct> selectPivots();

    /**
     * @brief Get distance from sparse matrix object
     *
     * @param pivotIndex1
     * @param pivotIndex2
     * @return float const
     */
    float const getDistance(unsigned int const pivotIndex1, unsigned int const pivotIndex2) {
        return _sparseMatrix->getDistance(pivotIndex1, pivotIndex2);
    }

  private:
    unsigned int const _datasetSize = 0;
    std::vector<float> const _pivotRadiusVector = std::vector<float>{};
    int const _numberOfPivotLayers = 1;

    std::shared_ptr<SparseMatrix> const _sparseMatrix;
};

#endif  // ConstantRadiusPivotSelection_hpp