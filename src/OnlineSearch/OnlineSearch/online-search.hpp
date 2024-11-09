#ifndef Search_hpp
#define Search_hpp
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
// 2022-05-31
#include <omp.h>
#include <chrono>

#include "pivot-layer-hierarchy.hpp"
#include "query-struct-online-search.hpp"
#include "search-sparse-matrix.hpp"

/**
 * @brief Class used to perform online search: find RNG neighbors of queries
 *
 */
class OnlineSearch {
  public:
    OnlineSearch(std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);
    ~OnlineSearch(){};
    inline void setNumThreads(int numThreads) { _numThreads = std::min(numThreads, omp_get_max_threads()); }

    void search(float* const& testPointer, unsigned int const testsetSize, unsigned int const dimension, bool verbose,
                tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& searchNeighbors);

    void _search(int const pivotLayerIndex, unsigned int const queryIndex, float const queryRadius, tsl::sparse_set<unsigned int>& queryNeighbors);

    void findNeighborsBruteForce(QueryStructOnlineSearch& queryStruct);

    void stage0(QueryStructOnlineSearch& queryStruct);
    void stage1(QueryStructOnlineSearch& queryStruct);
    void stage2(QueryStructOnlineSearch& queryStruct);
    void stage3(QueryStructOnlineSearch& queryStruct);
    void stage4(QueryStructOnlineSearch& queryStruct);
    void stage5(QueryStructOnlineSearch& queryStruct);
    void stage6(QueryStructOnlineSearch& queryStruct);

    void stage3_multithreading(QueryStructOnlineSearch& queryStruct);
    void stage6_multithreading(QueryStructOnlineSearch& queryStruct);

    void stage6_bruteForce(QueryStructOnlineSearch& queryStruct);

    inline float const getDistance(unsigned int const index1, unsigned int const index2) { return _sparseMatrix->getDistance(index1, index2); }
    inline std::pair<bool, float> const getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) {
        return _sparseMatrix->getDistanceIfAvailable(index1, index2);
    }
    unsigned long long int get_distanceComputationCount() const { return _sparseMatrix->get_distanceComputationCount(); }

    inline float _luneRadius(float d12, float r1, float r2) { return _luneRadiusVector[_luneIndex](d12, r1, r2); }

  private:
    std::shared_ptr<PivotLayerHierarchy> const _pivotLayers;

    std::shared_ptr<SearchSparseMatrix> _sparseMatrix;
    unsigned int _testsetSize = 0;

    int const _numberOfPivotLayers = 1;
    unsigned int const _datasetSize = 0;
    int _numThreads = 1;

    std::string _luneType = "Berk";
    int _luneIndex = 0;  // Berk default
    inline void _setLuneType(std::string luneType) {
        if (_luneType == "Cole") {
            _luneIndex = 1;
        }
    };
    inline static float _berkLune(float d12, float r1, float r2) { return d12 - 2 * r1 - r2; }
    inline static float _coleLune(float d12, float r1, float r2) { return d12 - r2; }
    std::vector<std::function<float(float, float, float)>> _luneRadiusVector =
        std::vector<std::function<float(float, float, float)>>{_berkLune, _coleLune};

  public:
    /**
     * @brief initialize metrics for offline index construction
     *
     */
    inline void initializeMetrics() {
        _distanceCount_layers.resize(_numberOfPivotLayers);
        _time_layers.resize(_numberOfPivotLayers);
        _distanceCount_stages.resize(7);
        _time_stages.resize(7);
    }

    // number of pivots
    std::vector<unsigned int> _numberOfPivotsPerLayer = std::vector<unsigned int>{};

    // distance computations
    unsigned long long int _distanceCount_build = 0;
    std::vector<unsigned long long int> _distanceCount_layers = std::vector<unsigned long long int>{};
    std::vector<unsigned long long int> _distanceCount_stages = std::vector<unsigned long long int>{};
    float _distanceCount_averageQuery = 0;

    // timing
    std::chrono::duration<double> _time_build;
    std::vector<std::chrono::duration<double>> _time_layers = std::vector<std::chrono::duration<double>>{};
    std::vector<std::chrono::duration<double>> _time_stages = std::vector<std::chrono::duration<double>>{};
    std::chrono::duration<double> _time_averageQuery;

    void printMetrics();
    void printSet(tsl::sparse_set<unsigned int> const& set);
};

#endif  // Search_hpp