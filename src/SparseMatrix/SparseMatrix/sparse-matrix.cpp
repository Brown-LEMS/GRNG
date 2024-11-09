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
#include "sparse-matrix.hpp"
#include <cmath>

SparseMatrix::SparseMatrix(float* const& dataPointer, unsigned int const datasetSize, unsigned int const dimension)
    : _dataPointer(dataPointer), _datasetSize(datasetSize), _dimension(dimension) {
    // initialize the sparse matrix structure
    _distanceMatrix.clear();
    _distanceMatrix.resize(_datasetSize);
    _distanceMatrix.shrink_to_fit();
}

SparseMatrix::SparseMatrix(float* const& dataPointer, unsigned int const datasetSize, unsigned int const dimension, bool cacheAll)
    : _dataPointer(dataPointer), _datasetSize(datasetSize), _dimension(dimension), _cacheAll(cacheAll) {
    // initialize the sparse matrix structure
    _distanceMatrix.clear();
    _distanceMatrix.resize(_datasetSize);
    _distanceMatrix.shrink_to_fit();
}

/**
 * @brief retrieve didstance between two indices, computes and stores if necessary. Modified to account for NOT caching all distances. If
 * _cacheALL==true, then all distances are cached, performs as normal. Otherwise, we use a special query distance to store all query distances. Sparse
 * matrix must manually be updated via query matrix.
 *
 * @param index1
 * @param index2
 * @return float const
 */
float const SparseMatrix::getDistance(unsigned int const index1, unsigned int const index2) {
    std::pair<bool, float> resultPair = getDistanceIfAvailable(index1, index2);

    if (!resultPair.first) {
        resultPair.second = _computeDistance(index1, index2);

        if (_cacheAll) {
            _insertDistanceIntoMatrix(index1, index2, resultPair.second);
        } else {
            if (index1 == _queryIndex) {
                _insertDistanceIntoQueryMatrix(index2, resultPair.second);
            } else if (index2 == _queryIndex) {
                _insertDistanceIntoQueryMatrix(index1, resultPair.second);
            }
        }
    }

    return resultPair.second;
}

/**
 * @brief compute distance between two indices
 * 
 * @param index1 
 * @param index2 
 * @return float const 
 */
float const SparseMatrix::_computeDistance(unsigned int const index1, unsigned int const index2) {
    _distanceComputationCount++;

    float distance = 0;
    for (unsigned int d = 0; d < _dimension; d++) {
        float difference = (_dataPointer[index1 * _dimension + d]) - (_dataPointer[index2 * _dimension + d]);
        distance += difference * difference;
    }

    return std::sqrt(distance);
}


/**
 * @brief retrieves distance iff already computed. if caching all, then check main matrix. otherwise, check query matrix for query distnaces and
 * otherwise main
 *
 * @param index1
 * @param index2
 * @return std::pair<bool, float> const
 */
std::pair<bool, float> const SparseMatrix::getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) const {
    std::pair<bool, float> resultPair = std::pair<bool, float>{false, INFINITY};

    if (_cacheAll) {
        resultPair = getDistanceFromMatrixIfAvailable(index1, index2);
    } else {
        if (index1 == _queryIndex) {
            resultPair = getDistanceFromQueryMatrixIfAvailable(index2);
        } else if (index2 == _queryIndex) {
            resultPair = getDistanceFromQueryMatrixIfAvailable(index1);
        } else {
            resultPair = getDistanceFromMatrixIfAvailable(index1, index2);
        }
    }

    return resultPair;
}

/**
 * @brief Check for and retrieve if distance is available in the main sparse matrix.
 *
 * @param index1
 * @param index2
 * @return std::pair<bool, float> const
 */
std::pair<bool, float> const SparseMatrix::getDistanceFromMatrixIfAvailable(unsigned int const index1, unsigned int const index2) const {
    std::pair<bool, float> resultPair = std::pair<bool, float>{false, INFINITY};

    ska::bytell_hash_map<unsigned int, float>::const_iterator index_position;
    if (index1 == index2) {
        resultPair.first = true;
        resultPair.second = 0.0f;
    } else {  // distance would be stored in the larger value
        if (index1 > index2) {
            index_position = _distanceMatrix[index1].find(index2);
            if (index_position != _distanceMatrix[index1].end()) {  // in
                resultPair.first = true;
                resultPair.second = index_position->second;
            }
        } else {
            index_position = _distanceMatrix[index2].find(index1);
            if (index_position != _distanceMatrix[index2].end()) {  // in
                resultPair.first = true;
                resultPair.second = index_position->second;
            }
        }
    }

    return resultPair;
}

/**
 * @brief Get distance from the query Sparse matrix if its available. simple map
 *
 * @param index
 * @return std::pair<bool, float> const
 */
std::pair<bool, float> const SparseMatrix::getDistanceFromQueryMatrixIfAvailable(unsigned int const index) const {
    std::pair<bool, float> resultPair = std::pair<bool, float>{false, INFINITY};

    tsl::sparse_map<unsigned int, float>::const_iterator index_position;
    if (index == _queryIndex) {
        resultPair.first = true;
        resultPair.second = 0.0f;
    } else {
        index_position = _queryDistanceMatrix.find(index);
        if (index_position != _queryDistanceMatrix.end()) {  // in
            resultPair.first = true;
            resultPair.second = index_position->second;
        }
    }

    return resultPair;
}


void SparseMatrix::_insertDistanceIntoMatrix(unsigned int const index1, unsigned int const index2, float const distance) {
    if (index1 > index2) {
        _distanceMatrix[index1].insert({index2, distance});
    } else if (index1 < index2) {
        _distanceMatrix[index2].insert({index1, distance});
    }
    _cachedDistanceCount++;
}

/**
 * @brief insert distance to index in query's sparse matrix
 *
 * @param index
 * @param distance
 */
void SparseMatrix::_insertDistanceIntoQueryMatrix(unsigned int const index, float const distance) { _queryDistanceMatrix.insert({index, distance}); }



/**
 * @brief Update SparseMatrix with specified indices
 *
 * @param indices
 */
void SparseMatrix::updateQueryMatrix(tsl::sparse_set<unsigned int> const& indices) {
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = indices.begin(); it1 != indices.end(); it1++) {
        unsigned int index = (*it1);
        float const distance = getDistance(_queryIndex,index);
        _insertDistanceIntoMatrix(_queryIndex, index, distance);
    }
}


//==================================================================================
//
//
//
//
//
//
//
//
//          Local Matrix for Multithreading
//
//
//
//
//
//
//
//
//==================================================================================



/**
 * @brief Return a local SparseMatrix copy for use in parallel processing
 *
 * @return SparseMatrixLocal
 */
SparseMatrixLocal SparseMatrix::getLocalCopy() const { return SparseMatrixLocal(_dataPointer, _datasetSize, _dimension); }

/**
 * @brief Get distance referencing SparseMatrix and given Local one
 *
 * @param index1
 * @param index2
 * @param localSparseMatrix
 * @return float const
 */
float const SparseMatrix::getDistance_LocalAndGlobal(unsigned int const index1, unsigned int const index2,
                                                     SparseMatrixLocal& localSparseMatrix) const {
    std::pair<bool, float> resultPair;
    std::pair<bool, float> const resultPair_global = getDistanceIfAvailable(index1, index2);
    if (resultPair_global.first) {
        return resultPair_global.second;
    } else {
        float const distance = localSparseMatrix.getDistance(index1, index2);
        return distance;
    }

    return resultPair.second;
}

/**
 * @brief get distance if computed, return bool==false if not. check global and local matrices
 *
 * @param index1
 * @param index2
 * @param localSparseMatrix
 * @return std::pair<bool, float> const
 */
std::pair<bool, float> const SparseMatrix::getDistanceIfAvailable_LocalAndGlobal(unsigned int const index1, unsigned int const index2,
                                                                                 SparseMatrixLocal& localSparseMatrix) const {
    std::pair<bool, float> resultPair = {0, 0.0f};
    std::pair<bool, float> const resultPair_global = getDistanceIfAvailable(index1, index2);

    if (resultPair_global.first) {
        return resultPair_global;
    } else {
        std::pair<bool, float> const resultPair_local = localSparseMatrix.getDistanceIfAvailable(index1, index2);
        if (resultPair_local.first) {
            return resultPair_local;
        }
    }

    return resultPair;
}

/**
 * @brief Add distances from local sparse matrix into global
 *
 * @param localSparseMatrix
 */
void SparseMatrix::updateSparseMatrixWithLocal(SparseMatrixLocal& localSparseMatrix) {
    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>> const& localMatrix = localSparseMatrix.getSparseMatrix();

    // add each distance computation into the global sparse matrix
    tsl::sparse_map<unsigned int, tsl::sparse_map<unsigned int, float>>::const_iterator it1;
    for (it1 = localMatrix.begin(); it1 != localMatrix.end(); it1++) {
        unsigned int const index1 = (*it1).first;
        tsl::sparse_map<unsigned int, float> const& index1LocalDistancesMap = (*it1).second;

        tsl::sparse_map<unsigned int, float>::const_iterator it2;
        for (it2 = index1LocalDistancesMap.begin(); it2 != index1LocalDistancesMap.end(); it2++) {
            unsigned int const index2 = (*it2).first;
            float const distance = (*it2).second;

            // store all distances
            if (_cacheAll) {
                _insertDistanceIntoMatrix(index1, index2, distance);

                // only keep query distances, store in query matrix
            } else {
                if (index1 == _queryIndex) {
                    _insertDistanceIntoQueryMatrix(index2, distance);
                } else if (index2 == _queryIndex) {
                    _insertDistanceIntoQueryMatrix(index1, distance);
                }
            }
        }
    }

    _distanceComputationCount += localSparseMatrix.getDistanceComputationCount();
    return;
}
