#ifndef IncrementalBuild_hpp
#define IncrementalBuild_hpp
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
// 2022-05-29
#include <memory.h>
#include <omp.h>

#include <chrono>

#include "pivot-layer-hierarchy.hpp"
#include "query-struct-incremental-build.hpp"

/**
 * @brief Class used to perform true incremental build
 *
 */
class IncrementalBuild {
  public:
    IncrementalBuild(std::shared_ptr<PivotLayerHierarchy> const& pivotLayers);
    ~IncrementalBuild(){};
    inline void setNumThreads(int numThreads) { _numThreads = std::min(numThreads, omp_get_max_threads()); }

    void build(bool verbose);
    void incrementalUpdate(unsigned int const queryIndex);

    void findNeighborsBruteForce(QueryStructIncrementalBuild& queryStruct);
    void findInterferenceBruteForce(QueryStructIncrementalBuild& queryStruct);

    void stage0(QueryStructIncrementalBuild& queryStruct);
    void stage1(QueryStructIncrementalBuild& queryStruct);
    void stage2(QueryStructIncrementalBuild& queryStruct);
    void stage3(QueryStructIncrementalBuild& queryStruct);
    void stage4(QueryStructIncrementalBuild& queryStruct);
    void stage5(QueryStructIncrementalBuild& queryStruct);
    void stage6(QueryStructIncrementalBuild& queryStruct);
    void stage7(QueryStructIncrementalBuild& queryStruct);
    void updateDataStructures(QueryStructIncrementalBuild& queryStruct);

    void stage8(QueryStructIncrementalBuild& queryStruct);
    void stage9(QueryStructIncrementalBuild& queryStruct);

    void stage3_multithreading(QueryStructIncrementalBuild& queryStruct);
    void stage6_multithreading(QueryStructIncrementalBuild& queryStruct);
    void stage7_multithreading(QueryStructIncrementalBuild& queryStruct);
    void stage8_multithreading(QueryStructIncrementalBuild& queryStruct);
    void stage9_multithreading(QueryStructIncrementalBuild& queryStruct);

    void stage6_bruteForce(QueryStructIncrementalBuild& queryStruct);
    void stage7_bruteForce(QueryStructIncrementalBuild& queryStruct);
    void stage8_bruteForce(QueryStructIncrementalBuild& queryStruct);
    void stage9_bruteForce(QueryStructIncrementalBuild& queryStruct);

    inline float const getDistance(unsigned int const index1, unsigned int const index2) { return _sparseMatrix->getDistance(index1, index2); }
    inline std::pair<bool, float> const getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) {
        return _sparseMatrix->getDistanceIfAvailable(index1, index2);
    }
    unsigned long long int get_distanceComputationCount() const { return _sparseMatrix->get_distanceComputationCount(); }

    inline float _luneRadius(float d12, float r1, float r2) { return _luneRadiusVector[_luneIndex](d12, r1, r2); }

  private:
    std::shared_ptr<PivotLayerHierarchy> const _pivotLayers;
    std::shared_ptr<SparseMatrix> const _sparseMatrix;
    bool _cacheAll = true;

    std::vector<float> _pivotRadiusVector = std::vector<float>{0};
    int const _numberOfPivotLayers = 1;
    unsigned int const _datasetSize = 0;
    std::string _luneType = "Berk";
    int _numThreads = 1;

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
        _distanceCount_stages.resize(11);
        _time_stages.resize(11);
    }

    // number of pivots
    std::vector<unsigned int> _numberOfPivotsPerLayer = std::vector<unsigned int>{};
    std::vector<unsigned int> _numberLinksPerLayer = std::vector<unsigned int>{};
    std::vector<float> _averageDegreePerLayer = std::vector<float>{};
    std::vector<unsigned int> _minDegreePerLayer = std::vector<unsigned int>{};
    std::vector<unsigned int> _maxDegreePerLayer = std::vector<unsigned int>{};
    std::size_t _currentRSS = 0;

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
    inline void printSet(tsl::sparse_set<unsigned int> const& set) {
      printf("{,");
      tsl::sparse_set<unsigned int>::const_iterator it1;
      for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,",*it1);
      }
      printf("}");
    }
};

#endif  // IncrementalBuild_hpp