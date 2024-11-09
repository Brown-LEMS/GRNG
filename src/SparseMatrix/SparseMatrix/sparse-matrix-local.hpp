#ifndef SparseMatrixLocal_hpp
#define SparseMatrixLocal_hpp
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
// 05-27-2022
#include "tsl/sparse_map.h"

/**
 * @brief Essentially a lightweight copy of the SparseMatrix. Used in multithreading to avoid parallel data saving collisions.
 *
 */
class SparseMatrixLocal {
  public:
    SparseMatrixLocal() {};
    SparseMatrixLocal(float* const& dataPointer, unsigned int const datasetSize, unsigned int const dimension);
    ~SparseMatrixLocal() {};


    float const _computeDistance(unsigned int const index1, unsigned int const index2);
    float const getDistance(unsigned int const index1, unsigned int const index2);
    std::pair<bool, float> const getDistanceIfAvailable(unsigned int const index1, unsigned int const index2);

    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>> const& getSparseMatrix() const {
        return _distanceMatrix;
    }
    inline unsigned long long int const getDistanceComputationCount() const {
        return _distanceComputationCount;
    };

  private:
    // data information
    float* const _dataPointer = NULL;
    unsigned int const _datasetSize = 0;
    unsigned int const _dimension = 0;

    // distances
    void _insertDistanceToMatrix(unsigned int const index1, unsigned int const index2, float const distance);
    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>> _distanceMatrix{};
    unsigned long long int _distanceComputationCount = 0;
};

#endif  // SparseMatrixLocal_hpp