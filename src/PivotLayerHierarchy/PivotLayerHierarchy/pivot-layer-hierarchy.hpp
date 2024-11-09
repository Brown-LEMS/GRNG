#ifndef PivotLayerHierarchy_hpp
#define PivotLayerHierarchy_hpp
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
// 01-26-2022

#include <vector>

#include "pivot-layer.hpp"

/**
 * @brief Class to handle all datastructures for the hierarchy
 *
 */
class PivotLayerHierarchy {
  public:
    PivotLayerHierarchy(){};
    PivotLayerHierarchy(std::vector<float> const pivotRadiusVector, std::string luneType, float* dataPointer, unsigned int const datasetSize,
                        unsigned int const dimension);
    PivotLayerHierarchy(std::vector<float> const pivotRadiusVector, std::string luneType, float* dataPointer, unsigned int const datasetSize,
                        unsigned int const dimension, bool cacheAll);
    PivotLayerHierarchy(std::vector<float> const pivotRadiusVector, std::string luneType, std::shared_ptr<SparseMatrix> sparseMatrix);
    PivotLayerHierarchy(PivotLayerHierarchy const& ref);
    ~PivotLayerHierarchy(){};

    void initialize_hierarchy();

    /**
     * @brief returns the specific pivotLayer object from our vector
     *
     * @param layerIndex
     * @return PivotLayer&
     */
    inline PivotLayer& operator[](int const layerIndex) { return *_pivotLayers[layerIndex]; }

    inline float const getDistance(unsigned int const index1, unsigned int const index2) { return _sparseMatrix->getDistance(index1, index2); }
    inline std::pair<bool, float> const getDistanceIfAvailable(unsigned int const index1, unsigned int const index2) {
        return _sparseMatrix->getDistanceIfAvailable(index1, index2);
    }
    unsigned long long int get_distanceComputationCount() const { return _sparseMatrix->get_distanceComputationCount(); }

    /**
     * @brief get the link distance of by perspective of p1
     *
     * @param d12 distance between p1,p2
     * @param r1 radius of p1's domain
     * @param r2 radius of p2's domain
     * @return float
     */
    inline float _luneRadius(float d12, float r1, float r2) { return _luneRadiusVector[_luneIndex](d12, r1, r2); }

  private:
    std::vector<float> _pivotRadiusVector = std::vector<float>{0.0f};
    int _numberOfPivotLayers = 1;
    std::string _luneType = "Berk";
    bool _verbose = false;

    // vector of functions for lune type
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

    std::vector<std::shared_ptr<PivotLayer>> _pivotLayers = std::vector<std::shared_ptr<PivotLayer>>{};
    std::shared_ptr<SparseMatrix> _sparseMatrix = std::make_shared<SparseMatrix>();

  public:
    bool coarseNeighborsAvailable = false;

    std::vector<float> const& get_pivotRadiusVector() const { return _pivotRadiusVector; }
    int const get_numberOfPivotLayers() const { return _numberOfPivotLayers; }
    std::string const& get_luneType() const { return _luneType; }
    std::shared_ptr<SparseMatrix> const& get_sparseMatrix() const { return _sparseMatrix; };
    unsigned int const get_datasetSize() const { return _sparseMatrix->get_datasetSize(); }
    unsigned int const get_dimension() const { return _sparseMatrix->get_dimension(); }
    inline std::vector<std::shared_ptr<PivotLayer>> get_pivotLayers() const& { return _pivotLayers; }
    inline std::vector<unsigned int> get_numberOfPivotsPerLayer() const {
        std::vector<unsigned int> numberOfPivotsPerLayer(_numberOfPivotLayers, 0);
        for (int i = 0; i < _numberOfPivotLayers; i++) numberOfPivotsPerLayer[i] = (*_pivotLayers[i]).get_numberOfPivots();
        return numberOfPivotsPerLayer;
    }
    void getDegreeStats(std::vector<unsigned int>& numLinks, std::vector<float>& aveDegree, std::vector<unsigned int>& minDegree,
                         std::vector<unsigned int>& maxDegree);
    void getRNG(tsl::sparse_map<unsigned int,tsl::sparse_set<unsigned int>>& neighbors) {
        neighbors = *(*_pivotLayers[_numberOfPivotLayers-1]).get_pivotLayerNeighbors_ptr();
    }

    inline void set_pivotRadiusVector(std::vector<float> pivotRadiusVector) {
        _pivotRadiusVector = pivotRadiusVector;
        _numberOfPivotLayers = pivotRadiusVector.size();
    }
    inline void set_luneType(std::string luneType) { _luneType = luneType; }
    inline void set_sparseMatrix(std::shared_ptr<SparseMatrix> sparseMatrix) { _sparseMatrix = sparseMatrix; }

    void setPivotIndicesAndRadii(int const layerIndex, tsl::sparse_set<unsigned int> const& pivotIndices,
                                 tsl::sparse_map<unsigned int, float> const& pivotRadii);
    void addPivotToLayer(int const layerIndex, unsigned int const pivotIndex, float const pivotRadius);

    void addPivotLayerNeighbors(int const layerIndex, unsigned int const pivotIndex1, tsl::sparse_set<unsigned int> const& pivotLayerNeighbors);
    void removePivotLayerNeighbors(int const layerIndex, tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& invalidatedPivotLinks);

    void updateAncestors(int const layerIndex, unsigned int const pivotIndex,
                         tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap);
    void updateDescendants(int const layerIndex, unsigned int const pivotIndex,
                           tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& descendantIndicesMap);
    void checkToUpdateDMax(int const layerIndex, unsigned int const pivotIndex, int const ancestorLayerIndex, unsigned int const ancestorIndex);

    void checkAncestorsToRecomputeUmaxDescendant(int const layerIndex, unsigned int const pivotIndex1, unsigned int const pivotIndex2);
    void updateUmaxDescendantMaxLinkDistance_(int const layerIndex, unsigned int const pivotIndex);
    void _computeUmaxDescendantMaxLinkDistance(int const layerIndex, unsigned int const ancestorIndex, int const finerLayerIndex);
    void checkToRecomputeUmaxMaximumOfDescendant(int const layerIndex, unsigned int const pivotIndex1, unsigned int const pivotIndex2);

    void addCoarsePivotLayerNeighbors(int const layerIndex, unsigned int const pivotIndex,
                                      tsl::sparse_set<unsigned int> const& coarsePivotLayerNeighbors);

    void addFinerCoarsePivotLayerNeighbors(int const layerIndex, unsigned int const pivotIndex,
                                           tsl::sparse_set<unsigned int> const& finerPivotLayerNeighbors);
    void removeFinerCoarsePivotLayerNeighbors(int const layerIndex,
                                              tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& invalidatedPivotLinks);
    
    void _updateVmaxMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex, unsigned int const coarsePivotIndex);
    void _checkToRecomputeVmaxMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex, unsigned int const coarseIndex);
    void _computeVmaxMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex);

    void checkAncestorsToRecomputeVmaxDescendantCoarse(int const layerIndex, unsigned int const finePivot, unsigned int const coarsePivot);
    void _updateVmaxDescendantMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex);
    void _computeVmaxDescendantMaxCoarseLinkDistance(int const layerIndex, unsigned int const pivotIndex, int const finerLayerIndex);
    void _checkToRecomputeVmaxMaximumOfDescendantCoarse(int const layerIndex, unsigned int const fineIndex, unsigned int const coarseIndex);
};

#endif  // PivotLayerHierarchy_hpp