#ifndef hRNGResults_hpp
#define hRNGResults_hpp
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
// 2022-06-04
#include <limits>
#include <vector>

/**
 * @brief Struct to hold all information after running hRNG as function evalutation. Used for saving progress of optimization.
 *
 */
struct hRNGResults {
    hRNGResults(){};
    hRNGResults(const hRNGResults& results) {
        _buildCount = results._buildCount;
        _pivotRadiusVector = results._pivotRadiusVector;
        _effectivePivotRadiusVector = results._effectivePivotRadiusVector;
        _numberOfPivotsVector = results._numberOfPivotsVector;
        _searchDistances = results._searchDistances;
        _searchTime = results._searchTime;
        _buildDistances = results._buildDistances;
        _buildTime = results._buildTime;
        _memoryUsage = results._memoryUsage;
        _linkCountVector = results._linkCountVector;
        _averageDegreeVector = results._averageDegreeVector;
        _minDegreeVector = results._minDegreeVector;
        _maxDegreeVector = results._maxDegreeVector;
    };
    ~hRNGResults(){};

    // results
    unsigned int _buildCount = 0;
    std::vector<float> _pivotRadiusVector{};
    std::vector<float> _effectivePivotRadiusVector{};
    std::vector<unsigned int> _numberOfPivotsVector{};
    float _searchDistances = std::numeric_limits<float>::infinity();
    float _searchTime = 0;
    unsigned long long int _buildDistances = 0;
    float _buildTime = 0;
    float _averageDistances = std::numeric_limits<float>::infinity();
    float _averageTime = 0;
    std::vector<unsigned int> _linkCountVector{};
    std::vector<float> _averageDegreeVector{};
    std::vector<unsigned int> _minDegreeVector{};
    std::vector<unsigned int> _maxDegreeVector{};
    float _memoryUsage = 0;
};

#endif  // hRNGResults_hpp