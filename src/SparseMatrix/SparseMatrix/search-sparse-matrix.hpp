#ifndef SearchSparseMatrix_hpp
#define SearchSparseMatrix_hpp
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
// 2022-05-31
#include "search-sparse-matrix-local.hpp"
#include "sparse-matrix.hpp"

class SearchSparseMatrix {
  public:
    // constructor/destructor
    SearchSparseMatrix(){};
    SearchSparseMatrix(float* const& testPointer, unsigned int const testSetSize, unsigned int const dimension,
                       std::shared_ptr<SparseMatrix> const& sparseMatrix);
    ~SearchSparseMatrix(){};

    float const getDistance(unsigned int const index1, unsigned int const index2);
    std::pair<bool, float> const getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) const;
    float const _computeDistance(unsigned int const index1, unsigned int const index2);
    void _insertDistanceIntoMatrix(unsigned int const index1, unsigned int const index2, float const distance);

    SearchSparseMatrixLocal getLocalCopy() const;
    float const getDistance_LocalAndGlobal(unsigned int const index1, unsigned int const index2, SearchSparseMatrixLocal& localSparseMatrix) const;
    std::pair<bool, float> const getDistanceIfAvailable_LocalAndGlobal(unsigned int const index1, unsigned int const index2,
                                                                       SearchSparseMatrixLocal& localSparseMatrix) const;
    void updateSparseMatrixWithLocal(SearchSparseMatrixLocal& localSparseMatrix);

    inline unsigned long long int const get_distanceComputationCount() const { return _distanceComputationCount; }
    unsigned long long int _searchQueryComputationCount = 0;

  private:
    unsigned int _dimension = 0;

    // original dataset
    float* const& _dataPointer = NULL;
    unsigned int const _datasetSize = 0;

    // new search set
    float* _testPointer;
    unsigned int _testsetSize = 0;

    // distance storage
    std::shared_ptr<SparseMatrix> const _sparseMatrix = std::make_shared<SparseMatrix>();
    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>> _distanceMatrix{};

    // metrics
    unsigned long long int _distanceComputationCount = 0;
};

#endif  // SearchSparseMatrix_hpp