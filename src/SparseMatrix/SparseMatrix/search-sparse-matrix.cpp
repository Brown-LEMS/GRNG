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
#include "search-sparse-matrix.hpp"

SearchSparseMatrix::SearchSparseMatrix(float* const& testPointer, unsigned int const testsetSize, unsigned int const dimension,
                                       std::shared_ptr<SparseMatrix> const& sparseMatrix)
    : _testPointer(testPointer),
      _testsetSize(testsetSize),
      _dimension(dimension),
      _sparseMatrix(sparseMatrix),
      _dataPointer(sparseMatrix->_dataPointer),
      _datasetSize(sparseMatrix->_datasetSize) {
    
    sparseMatrix->setCacheAll(true);
    
    // prep the search distance matrix
    _distanceMatrix.rehash(_testsetSize);  // prep size of search matrix
    for (int i = _datasetSize; i < (int)(_datasetSize + _testsetSize); i++) {
        _distanceMatrix[i] = tsl::sparse_map<unsigned int, float>{};
    }
    return;
}

std::pair<bool, float> const SearchSparseMatrix::getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) const {
    std::pair<bool, float> resultPair;
    tsl::sparse_map<unsigned int, float>::const_iterator index_position;

    resultPair.first = false;
    resultPair.second = std::numeric_limits<float>::infinity();

    if (index1 == index2) {
        resultPair.first = true;
        resultPair.second = 0;
    } else {
        if (index1 > index2) {
            if (index1 >= _datasetSize) {
                index_position = _distanceMatrix.at(index1).find(index2);
                if (index_position != _distanceMatrix.at(index1).end())  // in
                {
                    resultPair.first = true;
                    resultPair.second = index_position->second;
                }
            } else {
                resultPair = _sparseMatrix->getDistanceIfAvailable(index1, index2);
            }
        } else {
            if (index2 >= _datasetSize) {
                index_position = _distanceMatrix.at(index2).find(index1);
                if (index_position != _distanceMatrix.at(index2).end())  // in
                {
                    resultPair.first = true;
                    resultPair.second = index_position->second;
                }
            } else {
                resultPair = _sparseMatrix->getDistanceIfAvailable(index2, index1);
            }
        }
    }

    return resultPair;
}

/**
 * @brief get distance between index1,index2
 *
 * @param index1
 * @param index2
 * @return float const
 */
float const SearchSparseMatrix::getDistance(unsigned int const index1, unsigned int const index2) {
    std::pair<bool, float> resultPair = getDistanceIfAvailable(index1, index2);

    if (resultPair.first == false) {
        resultPair.second = _computeDistance(index1, index2);
        _insertDistanceIntoMatrix(index1, index2, resultPair.second);
    }

    return resultPair.second;
}

/**
 * @brief insert distance into Search sparseMatrix if involves query, otherwise main sparseMatrix
 *
 * @param index1
 * @param index2
 * @param distance
 */
void SearchSparseMatrix::_insertDistanceIntoMatrix(unsigned int const index1, unsigned int const index2, float const distance) {
    if (index1 >= _datasetSize || index2 >= _datasetSize)  // make sure we are within bounds
    {
        if (index1 > index2) {
            _distanceMatrix[index1].insert({index2, distance});
        } else if (index1 < index2) {
            _distanceMatrix[index2].insert({index1, distance});
        }
    } else {
        _sparseMatrix->_insertDistanceIntoMatrix(index1, index2, distance);
    }

    return;
}

/**
 * @brief compute distance between any two indices. considers dataset indices and queries
 *
 * @param index1
 * @param index2
 * @return float const
 */
float const SearchSparseMatrix::_computeDistance(unsigned int const index1, unsigned int const index2) {
    _distanceComputationCount++;
    float p1[_dimension], p2[_dimension];

    // get p1 where index1 in dataset or test set
    if (index1 < _datasetSize) {
        for (unsigned int d = 0; d < _dimension; d++) {
            p1[d] = _dataPointer[index1 * _dimension + d];
        }
    } else {
        _searchQueryComputationCount++;
        unsigned int const testIndex1 = index1 - _datasetSize;
        for (unsigned int d = 0; d < _dimension; d++) {
            p1[d] = _testPointer[testIndex1 * _dimension + d];
        }
    }

    // get p2 where index2 in dataset or test set
    if (index2 < _datasetSize) {
        for (unsigned int d = 0; d < _dimension; d++) {
            p2[d] = _dataPointer[index2 * _dimension + d];
        }
    } else {
        _searchQueryComputationCount++;
        unsigned int const testIndex2 = index2 - _datasetSize;
        for (unsigned int d = 0; d < _dimension; d++) {
            p2[d] = _testPointer[testIndex2 * _dimension + d];
        }
    }

    // compute euclidean distance
    float distance = 0;
    for (unsigned int d = 0; d < _dimension; d++) {
        float difference = p1[d] - p2[d];
        distance += difference * difference;
    }

    return std::sqrt(distance);
}

/**
 *
 * Local Search Sparse Matrix functionality for Multithreading
 *
 */

SearchSparseMatrixLocal SearchSparseMatrix::getLocalCopy() const {
    return SearchSparseMatrixLocal(_dataPointer, _datasetSize, _testPointer, _testsetSize, _dimension);
}

/**
 * @brief Get distance between two indices. Check global sparseMatrix first, then consult local to search and compute if necessary
 *
 * @param index1
 * @param index2
 * @param localSparseMatrix
 * @return float const
 */
float const SearchSparseMatrix::getDistance_LocalAndGlobal(unsigned int const index1, unsigned int const index2,
                                                           SearchSparseMatrixLocal& localSparseMatrix) const {
    std::pair<bool, float> resultPair;

    // check global matrix first
    std::pair<bool, float> const resultPair_global = getDistanceIfAvailable(index1, index2);
    if (resultPair_global.first) {
        return resultPair_global.second;
    } else {  // check local next, computing if necessary
        float const distance = localSparseMatrix.getDistance(index1, index2);
        return distance;
    }

    return resultPair.second;
}

/**
 * @brief Get distance only if stored, check global and then local sparse matrix
 *
 * @param index1
 * @param index2
 * @param localSparseMatrix
 * @return std::pair<bool, float> const
 */
std::pair<bool, float> const SearchSparseMatrix::getDistanceIfAvailable_LocalAndGlobal(unsigned int const index1, unsigned int const index2,
                                                                                       SearchSparseMatrixLocal& localSparseMatrix) const {
    std::pair<bool, float> resultPair = {0, 0.0f};

    // try getting from global
    std::pair<bool, float> const resultPair_global = getDistanceIfAvailable(index1, index2);
    if (resultPair_global.first) {
        return resultPair_global;
    } else {  // try getting from local
        std::pair<bool, float> const resultPair_local = localSparseMatrix.getDistanceIfAvailable(index1, index2);
        if (resultPair_local.first) {
            return resultPair_local;
        }
    }

    return resultPair;
}

/**
 * @brief update global sparse matrix with all distances stored within local. Done at conclusion of multithreading.
 *
 * @param localSparseMatrix
 */
void SearchSparseMatrix::updateSparseMatrixWithLocal(SearchSparseMatrixLocal& localSparseMatrix) {
    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>> const& localMatrix = localSparseMatrix.get_sparseMatrix();

    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>>::const_iterator it1;
    for (it1 = localMatrix.begin(); it1 != localMatrix.end(); it1++) {
        unsigned int const index1 = (*it1).first;
        tsl::sparse_map<unsigned int, float> const& index1LocalDistancesMap = (*it1).second;

        tsl::sparse_map<unsigned int, float>::const_iterator it2;
        for (it2 = index1LocalDistancesMap.begin(); it2 != index1LocalDistancesMap.end(); it2++) {
            unsigned int const index2 = (*it2).first;
            float const distance = (*it2).second;
            _insertDistanceIntoMatrix(index1, index2, distance);
        }
    }
    _distanceComputationCount += localSparseMatrix.get_distanceComputationCount();
    return;
}