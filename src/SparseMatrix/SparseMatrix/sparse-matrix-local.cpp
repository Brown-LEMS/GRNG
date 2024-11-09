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
#include "sparse-matrix-local.hpp"

SparseMatrixLocal::SparseMatrixLocal(float* const& dataPointer, unsigned int const datasetSize, unsigned int const dimension)
    : _dataPointer(dataPointer), _datasetSize(datasetSize), _dimension(dimension){};

/**
 * @brief Get distance between two indices. Checks if stored first. Computes and stores if not.
 *
 * @param index1
 * @param index2
 * @return float const
 */
float const SparseMatrixLocal::getDistance(unsigned int const index1, unsigned int const index2) {
    std::pair<bool, float> resultPair = getDistanceIfAvailable(index1, index2);

    if (resultPair.first == false) {
        resultPair.second = _computeDistance(index1, index2);
        _insertDistanceToMatrix(index1, index2, resultPair.second);
    }

    return resultPair.second;
}

/**
 * @brief Get distance between two indices if stored. Return bool==false otherwise. const eventually?
 *
 * @param index1
 * @param index2
 * @return std::pair<bool, float> const
 */
std::pair<bool, float> const SparseMatrixLocal::getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) {
    std::pair<bool, float> resultPair;
    tsl::sparse_map<unsigned int, float>::const_iterator index_position;

    resultPair.first = false;
    resultPair.second = std::numeric_limits<float>::infinity();

    if (index1 == index2) {
        resultPair.first = true;
        resultPair.second = 0;
    } else {
        if (index1 > index2) {
            index_position = _distanceMatrix[index1].find(index2);
            if (index_position != _distanceMatrix[index1].end())  // in
            {
                resultPair.first = true;
                resultPair.second = index_position->second;
            }
        } else {
            index_position = _distanceMatrix[index2].find(index1);
            if (index_position != _distanceMatrix[index2].end())  // in
            {
                resultPair.first = true;
                resultPair.second = index_position->second;
            }
        }
    }

    return resultPair;
}

float const SparseMatrixLocal::_computeDistance(unsigned int const index1, unsigned int const index2) {
    _distanceComputationCount++;

    float distance = 0;
    for (unsigned int d = 0; d < _dimension; d++) {
        float difference = (_dataPointer[index1 * _dimension + d]) - (_dataPointer[index2 * _dimension + d]);
        distance += difference * difference;
    }

    return std::sqrt(distance);
}

/**
 * @brief insert distance into own sparse matrix, not worrying about dataset vs. testset.
 * @note is this safe? should we initialize maps before inserting?
 *
 * @param index1
 * @param index2
 * @param distance
 */
void SparseMatrixLocal::_insertDistanceToMatrix(unsigned int const index1, unsigned int const index2, float const distance) {
    if (index1 > index2) {
        _distanceMatrix[index1].insert({index2, distance});
    } else if (index1 < index2) {
        _distanceMatrix[index2].insert({index1, distance});
    }
}