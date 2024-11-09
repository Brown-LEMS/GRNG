#ifndef HierarchicalRNG_hpp
#define HierarchicalRNG_hpp
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
// 01-24-2022

#include "pivot-layer-hierarchy.hpp"
#include "sparse-matrix.hpp"

/**
 * @brief Main Interface for Hierarchical RNG
 *
 */
class HierarchicalRNG {
  public:
    HierarchicalRNG(float* const dataPointer, unsigned int const datasetSize, unsigned int const dimension);
    HierarchicalRNG(float* const dataPointer, unsigned int const datasetSize, unsigned int const dimension, bool cacheAll);
    HierarchicalRNG(std::shared_ptr<SparseMatrix> sparseMatrix);
    HierarchicalRNG(std::shared_ptr<PivotLayerHierarchy> pivotLayers);
    HierarchicalRNG(){};
    ~HierarchicalRNG(){};

    std::vector<float> getEffectiveRadiusVector(std::vector<float> const vec);

    void offlineBuild(std::vector<float> const pivotRadiusVector, std::string luneType, int numThreads, bool verbose);
    void incrementalBuild(std::vector<float> const pivotRadiusVector, std::string luneType, int numThreads, bool verbose);
    void onlineSearch(float* const testPointer, unsigned int const testsetSize, unsigned int const dimension, int numThreads, bool verbose);

  private:
    /**
     * @brief Pointer to the all Hierarchical Data Structures
     *
     */
    std::shared_ptr<PivotLayerHierarchy> _pivotLayers;
    std::shared_ptr<SparseMatrix> _sparseMatrix;
};

#endif  // HierarchicalRNG_hpp