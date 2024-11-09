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

#include "hierarchical-rng.hpp"

#include "incremental-build.hpp"
#include "offline-build.hpp"
#include "online-search.hpp"

/**
 * @brief General constructor for new hierarchy
 *
 * @param dataPointer // float pointer as a flattened NxD array
 * @param datasetSize // dataset size N
 * @param dimension   // dimension D
 */
HierarchicalRNG::HierarchicalRNG(float* const dataPointer, unsigned int const datasetSize, unsigned int const dimension) {
    _sparseMatrix = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension);
}

HierarchicalRNG::HierarchicalRNG(float* const dataPointer, unsigned int const datasetSize, unsigned int const dimension, bool cacheAll) {
    _sparseMatrix = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension, cacheAll);
}

/**
 * @brief Constructor for new hierarchy, giving a saved SparseMatrix
 *
 * @param sparseMatrix
 */
HierarchicalRNG::HierarchicalRNG(std::shared_ptr<SparseMatrix> sparseMatrix) { _sparseMatrix = sparseMatrix; }

/**
 * @brief Returns new pivot radius vector accounting for relative pivot domain
 *
 * @param vec
 * @return std::vector<float>
 */
std::vector<float> HierarchicalRNG::getEffectiveRadiusVector(std::vector<float> const vec) {
    std::vector<float> new_vec = vec;
    for (int i = 0; i < (int)vec.size(); i++) {
        for (int j = i + 1; j < (int)vec.size(); j++) {
            new_vec[i] += (float)vec[j];
        }
    }
    return new_vec;
}

/**
 * @brief Perform Offline Index Construction via Incremental with Pivot Selection.
 *
 * @param pivotRadiusVector
 * @param luneType // type of lune for hierarchy. "Berk" gives exact RNG.
 * @param numThreads
 * @param verbose
 */
void HierarchicalRNG::offlineBuild(std::vector<float> const pivotRadiusVector, std::string luneType, int numThreads, bool verbose) {
    _pivotLayers = std::make_shared<PivotLayerHierarchy>(pivotRadiusVector, luneType, _sparseMatrix);

    OfflineBuild OB(_pivotLayers);
    OB.setNumThreads(numThreads);
    OB.build(verbose);
}

void HierarchicalRNG::incrementalBuild(std::vector<float> const pivotRadiusVector, std::string luneType, int numThreads, bool verbose) {
    _pivotLayers = std::make_shared<PivotLayerHierarchy>(pivotRadiusVector, luneType, _sparseMatrix);

    IncrementalBuild IB(_pivotLayers);
    IB.setNumThreads(numThreads);
    IB.build(verbose);
}

void HierarchicalRNG::onlineSearch(float* const testPointer, unsigned int const testsetSize, unsigned int const dimension, int numThreads,
                                   bool verbose) {
    OnlineSearch OS(_pivotLayers);
    OS.setNumThreads(numThreads);
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> searchNeighbors;
    OS.search(testPointer, testsetSize, dimension, verbose, searchNeighbors);
}