#ifndef SparseMatrix_h
#define SparseMatrix_h
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

#include <vector>

#include "bytell_hash_map.hpp"
#include "sparse-matrix-local.hpp"
#include "tsl/sparse_set.h"
#include "tsl/sparse_map.h"

/**
 * @brief sparse matrix to hold all pairwise distance computations computed. Stored as lower triangular matrix. Uses a vector of hashmaps for speed.
 *
 */
class SparseMatrix {
  public:
    SparseMatrix(){};
    SparseMatrix(float *const &dataPointer, unsigned int const datasetSize, unsigned int const dimension);
    SparseMatrix(float *const &dataPointer, unsigned int const datasetSize, unsigned int const dimension, bool cacheAll);
    ~SparseMatrix(){};

    float const getDistance(unsigned int const index1, unsigned int const index2);

    inline void setNewQuery(unsigned int const queryIndex) {
        _queryIndex = queryIndex;
        _queryDistanceMatrix.clear();
    }
    void updateQueryMatrix(tsl::sparse_set<unsigned int> const& indices);

    std::pair<bool, float> const getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) const;
    std::pair<bool, float> const getDistanceFromMatrixIfAvailable(unsigned int const index1, unsigned int const index2) const;
    std::pair<bool, float> const getDistanceFromQueryMatrixIfAvailable(unsigned int const index) const;

    void _insertDistanceIntoMatrix(unsigned int const index1, unsigned int const index2, float const distance);
    void _insertDistanceIntoQueryMatrix(unsigned int const index, float const distance);

    float const _computeDistance(unsigned int const index1, unsigned int const index2);

    SparseMatrixLocal getLocalCopy() const;
    float const getDistance_LocalAndGlobal(unsigned int const index1, unsigned int const index2, SparseMatrixLocal &localSparseMatrix) const;
    std::pair<bool, float> const getDistanceIfAvailable_LocalAndGlobal(unsigned int const index1, unsigned int const index2,
                                                                       SparseMatrixLocal &localSparseMatrix) const;
    void updateSparseMatrixWithLocal(SparseMatrixLocal &localSparseMatrix);

    float *const &get_dataPointer() { return _dataPointer; }
    unsigned int const get_datasetSize() const { return _datasetSize; }
    unsigned int const get_dimension() const { return _dimension; }
    std::vector<ska::bytell_hash_map<unsigned int, float>> const &get_distanceMatrix() const { return _distanceMatrix; }
    unsigned long long int const get_distanceComputationCount() { return _distanceComputationCount; }

    inline bool getCacheAllStatus() const {return _cacheAll; }
    inline void setCacheAll(bool cacheAll) {
      _cacheAll = cacheAll;
    }

    // void set_dataPointer(float *const &dataPointer) { _dataPointer = dataPointer; }
    // void set_datasetSize(unsigned int const &datasetSize) { _datasetSize = datasetSize; }
    // void set_dimension(unsigned int const &dimension) { _dimension = dimension; }
    // void set_distanceMatrix(std::vector<ska::bytell_hash_map<unsigned int, float>> &distanceMatrix) { _distanceMatrix = distanceMatrix; }

    float *_dataPointer = NULL;
    unsigned int _datasetSize = 0;

  private:
    unsigned int _dimension = 0;

    // for caching only relationship distances
    unsigned int _queryIndex = 0;
    bool _cacheAll = true;

    // distance storage
    std::vector<ska::bytell_hash_map<unsigned int, float>> _distanceMatrix{};
    tsl::sparse_map<unsigned int, float> _queryDistanceMatrix{};

    // statistics
    unsigned long long int _distanceComputationCount = 0;
    unsigned long long int _cachedDistanceCount = 0;
};

#endif  // SparseMatrix_h