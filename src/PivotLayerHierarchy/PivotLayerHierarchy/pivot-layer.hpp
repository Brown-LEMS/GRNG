#ifndef PivotLayer_hpp
#define PivotLayer_hpp
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

#include <memory.h>

#include "index-pair-distance-struct.hpp"
#include "sparse-matrix.hpp"
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"

/**
 * @brief Class holding all datastructures for an individual Pivot Layer
 *
 */
class PivotLayer {
  public:
    PivotLayer(){};
    /**
     * @brief Construct a new Pivot Layer object
     *
     * @param pivotLayerID layer id within the hierarchy
     * @param numberOfPivotLayers number of layers in the hierarchy
     * @param SparseMatrix shared pointer to SparseMatrix object
     * @param luneType the lune type to use
     * @param verbose
     */
    PivotLayer(int const pivotLayerID, int const numberOfPivotLayers, std::shared_ptr<SparseMatrix> const SparseMatrix, std::string luneType,
               bool verbose);
    PivotLayer(PivotLayer const& layer);
    ~PivotLayer(){};

    /**
     * @brief Add a single pivot to the layer
     *
     * @param pivotIndex
     * @param radius
     */
    void addPivotToLayer(unsigned int const pivotIndex, float const radius);

    /**
     * @brief Initialize all data structures for a single pivot
     *
     * @param pivotIndex
     */
    void initializePivotInLayer(unsigned int const pivotIndex);

    /**
     * @brief add link(pivotIndex1,pivotIndex2) to the PivotLayerGraph
     *
     * @param pivotIndex1
     * @param pivotIndex2
     */
    void addLinkToLayer(unsigned int const pivotIndex1, unsigned int const pivotIndex2);

    /**
     * @brief remove link(pivotIndex1,pivotIndex2) from the PivotLayerGraph
     *
     * @param pivotIndex1
     * @param pivotIndex2
     */
    void removeLinkFromLayer(unsigned int const pivotIndex1, unsigned int const pivotIndex2);

    /**
     * @brief add ancestorIndex as an ancestor of pivotIndex in layer ancestorLayerIndex
     *
     * @param pivotIndex
     * @param ancestorLayerIndex
     * @param ancestorIndex
     */
    inline void addAncestor(unsigned int const pivotIndex, int const ancestorLayerIndex, unsigned int const ancestorIndex) {
        (*_ancestorPivotIndices)[pivotIndex][ancestorLayerIndex].insert(ancestorIndex);
    }

    /**
     * @brief add descendantIndex as a descendant of pivotIndex in layer descendantLayerIndex
     *
     * @param pivotIndex
     * @param descendantLayerIndex
     * @param descendantIndex
     */
    inline void addDescendant(unsigned int const pivotIndex, int const descendantLayerIndex, unsigned int const descendantIndex) {
        (*_descendantPivotIndices)[pivotIndex][descendantLayerIndex].insert(descendantIndex);
    }

    void _updateUmaxMaxLinkDistance(unsigned int const pivotIndex1, unsigned int const pivotIndex2);
    void _computeUmaxMaxLinkDistance(unsigned int const pivotIndex1);

    // VMAX
    void _updateVmaxMaxCoarseLinkDistance(unsigned int const pivotIndex, unsigned int const coarsePivotIndex, float const coarsePivotRadius);
    void _computeVmaxMaxCoarseLinkDistance(unsigned int const pivotIndex);

    /**
     * @brief compute the pivotLayerGraph by brute force. Used for top layer initialization.
     *
     */
    void computeBruteForcePivotLayerGraph();

    /**
     * @brief retrieve distance from SparseMatrix
     *
     * @param index1
     * @param index2
     * @return float const
     */
    inline float const getDistance(unsigned int const index1, unsigned int const index2) { return _sparseMatrix->getDistance(index1, index2); }

    /**
     * @brief Get the distanceComputationCount from SparseMatrix
     *
     * @return unsigned long long int
     */
    unsigned long long int get_distanceComputationCount() const { return _sparseMatrix->get_distanceComputationCount(); }

    unsigned int const get_numberOfPivots() const { return _pivotIndices->size(); }

    /**
     * @brief return the distance regarding lune type
     *
     * @param d12
     * @param r1
     * @param r2
     * @return float
     */
    inline float _luneRadius(float d12, float r1, float r2) { return _luneRadiusVector[_luneIndex](d12, r1, r2); }

  private:
    // constructed vars
    int const _pivotLayerID = 0;
    int const _numberOfPivotLayers = 1;
    std::shared_ptr<SparseMatrix> _sparseMatrix = std::make_shared<SparseMatrix>();
    std::string _luneType = "Berk";
    bool _verbose = false;

    // vector of functions for lune type
    int _luneIndex = 0;  // Berk default
    inline static float _berkLune(float d12, float r1, float r2) { return d12 - 2 * r1 - r2; }
    inline static float _coleLune(float d12, float r1, float r2) { return d12 - r2; }
    std::vector<std::function<float(float, float, float)>> _luneRadiusVector =
        std::vector<std::function<float(float, float, float)>>{_berkLune, _coleLune};

    // DATA STRUCTURES (as shared pointers)
    std::shared_ptr<tsl::sparse_set<unsigned int>> _pivotIndices = std::make_shared<tsl::sparse_set<unsigned int>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, float>> _pivotRadii = std::make_shared<tsl::sparse_map<unsigned int, float>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> _pivotLayerNeighbors =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> _coarsePivotLayerNeighbors =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>> _ancestorPivotIndices =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>> _descendantPivotIndices =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>>> _dmaxMaxChildDistance =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>>>();
    std::shared_ptr<tsl::sparse_map<int, float>> _dmaxMaximumOfMaxChildDistance = std::make_shared<tsl::sparse_map<int, float>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>> _umaxMaxLinkDistance =
        std::make_shared<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>> _umaxDescendantMaxLinkDistance =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>>();
    std::shared_ptr<tsl::sparse_map<int, IndexPairDistanceStruct>> _umaxMaximumOfDescendantMaxLinkDistance =
        std::make_shared<tsl::sparse_map<int, IndexPairDistanceStruct>>();

    // VMAX
    std::shared_ptr<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>> _vmaxMaxCoarseLinkDistance =
        std::make_shared<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>> _vmaxDescendantMaxCoarseLinkDistance =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>>();
    std::shared_ptr<tsl::sparse_map<int, IndexPairDistanceStruct>> _vmaxMaximumOfDescendantMaxCoarseLinkDistance =
        std::make_shared<tsl::sparse_map<int, IndexPairDistanceStruct>>();

  public:
    // general get/set
    int const get_pivotLayerID() const { return _pivotLayerID; }
    int const get_numberOfPivotLayers() const { return _numberOfPivotLayers; }
    std::shared_ptr<SparseMatrix> get_sparseMatrix() const { return _sparseMatrix; }
    std::string const get_luneType() const { return _luneType; }

    // GETTING FUNCTIONS (for build/search)
    inline float const get_pivotRadius(unsigned int const pivotIndex) const { return (*_pivotRadii).at(pivotIndex); }
    inline tsl::sparse_set<unsigned int> const& get_pivotLayerNeighbors(unsigned int const pivotIndex) const {
        return (*_pivotLayerNeighbors).at(pivotIndex);
    }
    inline tsl::sparse_set<unsigned int> const& get_coarsePivotLayerNeighbors(unsigned int const pivotIndex) const {
        return (*_coarsePivotLayerNeighbors).at(pivotIndex);
    }
    inline tsl::sparse_set<unsigned int> const& get_ancestorPivotIndices(unsigned int const pivotIndex, int const coarserLayerIndex) const {
        return (*_ancestorPivotIndices).at(pivotIndex).at(coarserLayerIndex);
    }
    inline tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& get_ancestorPivotIndicesMap(unsigned int const pivotIndex) const {
        return (*_ancestorPivotIndices).at(pivotIndex);
    }
    inline tsl::sparse_set<unsigned int> const& get_descendantPivotIndices(unsigned int const pivotIndex, int const finerLayerIndex) const {
        return (*_descendantPivotIndices).at(pivotIndex).at(finerLayerIndex);
    }
    inline tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& get_descendantPivotIndicesMap(unsigned int const pivotIndex) const {
        return (*_descendantPivotIndices).at(pivotIndex);
    }
    inline float const& get_dmaxMaxChildDistance(unsigned int const pivotIndex, int const finerLayerIndex) const {
        return (*_dmaxMaxChildDistance).at(pivotIndex).at(finerLayerIndex);
    }
    inline float const& get_dmaxMaximumOfMaxChildDistance(int const finerLayerIndex) const {
        return (*_dmaxMaximumOfMaxChildDistance).at(finerLayerIndex);
    }
    inline IndexPairDistanceStruct const& get_umaxMaxLinkDistance(unsigned int const pivotIndex) const {
        return (*_umaxMaxLinkDistance).at(pivotIndex);
    }
    inline IndexPairDistanceStruct const& get_umaxDescendantMaxLinkDistance(unsigned int const pivotIndex, int const finerLayerIndex) const {
        return (*_umaxDescendantMaxLinkDistance).at(pivotIndex).at(finerLayerIndex);
    }
    inline IndexPairDistanceStruct const& get_umaxMaximumOfDescendantMaxLinkDistance(int const finerLayerIndex) const {
        return (*_umaxMaximumOfDescendantMaxLinkDistance).at(finerLayerIndex);
    }

    // VMAX
    inline IndexPairDistanceStruct const& get_vmaxMaxCoarseLinkDistance(unsigned int const pivotIndex) const {
        return (*_vmaxMaxCoarseLinkDistance).at(pivotIndex);
    }
    inline IndexPairDistanceStruct const& get_vmaxDescendantMaxCoarseLinkDistance(unsigned int const pivotIndex, int const finerLayerIndex) const {
        return (*_vmaxDescendantMaxCoarseLinkDistance).at(pivotIndex).at(finerLayerIndex);
    }
    inline IndexPairDistanceStruct const& get_vmaxMaximumOfDescendantMaxCoarseLinkDistance(int const finerLayerIndex) const {
        return (*_vmaxMaximumOfDescendantMaxCoarseLinkDistance).at(finerLayerIndex);
    }

    // SETTING FUNCTIONS (for build/search)
    inline void set_coarsePivotLayerNeighbors(unsigned int const pivotIndex, tsl::sparse_set<unsigned int> const& coarseNeighbors) {
        (*_coarsePivotLayerNeighbors)[pivotIndex] = coarseNeighbors;
    }
    void add_coarsePivotLayerNeighbor(unsigned int const pivotIndex, unsigned int const coarseNeighbor);
    void remove_coarsePivotLayerNeighbor(unsigned int const pivotIndex, unsigned int const coarseNeighbor); 
    inline void set_dmaxMaxChildDistance(unsigned int const pivotIndex, int const finerLayerIndex, float const distance) {
        (*_dmaxMaxChildDistance)[pivotIndex][finerLayerIndex] = distance;
    }
    inline void set_dmaxMaximumOfMaxChildDistance(int const finerLayerIndex, float const distance) {
        (*_dmaxMaximumOfMaxChildDistance)[finerLayerIndex] = distance;
    }
    inline void set_umaxDescendantMaxLinkDistance(unsigned int const pivotIndex, int const finerLayerIndex, IndexPairDistanceStruct const& umax) {
        (*_umaxDescendantMaxLinkDistance)[pivotIndex][finerLayerIndex] = umax;
    }
    inline void set_umaxMaximumOfDescendantMaxLinkDistance(int const finerLayerIndex, IndexPairDistanceStruct const& umax) {
        (*_umaxMaximumOfDescendantMaxLinkDistance)[finerLayerIndex] = umax;
    }

    // VMAX
    inline void set_vmaxMaxCoarseLinkDistance(unsigned int const pivotIndex, IndexPairDistanceStruct const& vmax) {
        (*_vmaxMaxCoarseLinkDistance)[pivotIndex] = vmax;
    }
    inline void set_vmaxDescendantMaxCoarseLinkDistance(unsigned int const pivotIndex, int const finerLayerIndex, IndexPairDistanceStruct const& vmax) {
        (*_vmaxDescendantMaxCoarseLinkDistance)[pivotIndex][finerLayerIndex] = vmax;
    }
    inline void set_vmaxMaximumOfDescendantMaxCoarseLinkDistance(int const finerLayerIndex, IndexPairDistanceStruct const& vmax) {
        (*_vmaxMaximumOfDescendantMaxCoarseLinkDistance)[finerLayerIndex] = vmax;
    }

    // GETTING FUNCTIONS (for shared ptrs, copy constructor)
    inline std::shared_ptr<tsl::sparse_set<unsigned int>> get_pivotIndices_ptr() const { return _pivotIndices; }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, float>> get_pivotRadii_ptr() const { return _pivotRadii; }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> get_pivotLayerNeighbors_ptr() const {
        return _pivotLayerNeighbors;
    }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> get_coarsePivotLayerNeighbors_ptr() const {
        return _coarsePivotLayerNeighbors;
    }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>> get_ancestorPivotIndices_ptr() const {
        return _ancestorPivotIndices;
    }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>> get_descendantPivotIndices_ptr()
        const {
        return _descendantPivotIndices;
    }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>>> get_dmaxMaxChildDistance_ptr() const {
        return _dmaxMaxChildDistance;
    }
    inline std::shared_ptr<tsl::sparse_map<int, float>> get_dmaxMaximumOfMaxChildDistance_ptr() const { return _dmaxMaximumOfMaxChildDistance; }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>> get_umaxMaxLinkDistance_ptr() const {
        return _umaxMaxLinkDistance;
    }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>> get_umaxDescendantMaxLinkDistance_ptr()
        const {
        return _umaxDescendantMaxLinkDistance;
    }
    inline std::shared_ptr<tsl::sparse_map<int, IndexPairDistanceStruct>> get_umaxMaximumOfDescendantMaxLinkDistance_ptr() const {
        return _umaxMaximumOfDescendantMaxLinkDistance;
    }

    // VMAX
    inline std::shared_ptr<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>> get_vmaxMaxCoarseLinkDistance_ptr() const {
        return _vmaxMaxCoarseLinkDistance;
    }
    inline std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>> get_vmaxDescendantMaxCoarseLinkDistance_ptr()
        const {
        return _vmaxDescendantMaxCoarseLinkDistance;
    }
    inline std::shared_ptr<tsl::sparse_map<int, IndexPairDistanceStruct>> get_vmaxMaximumOfDescendantMaxCoarseLinkDistance_ptr() const {
        return _vmaxMaximumOfDescendantMaxCoarseLinkDistance;
    }
    

    // SETTING FUNCTIONS (for deserialization)
    inline void set_pivotIndices_ptr(std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices) { _pivotIndices = pivotIndices; }
    inline void set_pivotLayerRadii_ptr(std::shared_ptr<tsl::sparse_map<unsigned int, float>> pivotRadii) { _pivotRadii = pivotRadii; }
    inline void set_pivotLayerNeighbors_ptr(std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> pivotLayerNeighbors) {
        _pivotLayerNeighbors = pivotLayerNeighbors;
    }
    inline void set_coarsePivotLayerNeighbors_ptr(
        std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> coarserPivotLayerNeighbors) {
        _coarsePivotLayerNeighbors = coarserPivotLayerNeighbors;
    }
    inline void set_ancestorPivotIndices_ptr(
        std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>> ancestorPivotIndices) {
        _ancestorPivotIndices = ancestorPivotIndices;
    }
    inline void set_descendantPivotIndices_ptr(
        std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>>> descendantPivotIndices) {
        _descendantPivotIndices = descendantPivotIndices;
    }
    inline void set_dmaxMaxChildDistance_ptr(std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>>> dmaxMaxChildDistance) {
        _dmaxMaxChildDistance = dmaxMaxChildDistance;
    }
    inline void set_dmaxMaximumOfMaxChildDistance_ptr(std::shared_ptr<tsl::sparse_map<int, float>> dmaxMaximumOfMaxChildDistance) {
        _dmaxMaximumOfMaxChildDistance = dmaxMaximumOfMaxChildDistance;
    }
    inline void set_umaxMaxLinkDistance_ptr(std::shared_ptr<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>> umaxMaxLinkDistance) {
        _umaxMaxLinkDistance = umaxMaxLinkDistance;
    }
    inline void set_umaxDescendantMaxLinkDistance_ptr(
        std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>> umaxDescendantMaxLinkDistance) {
        _umaxDescendantMaxLinkDistance = umaxDescendantMaxLinkDistance;
    }
    inline void set_umaxMaximumOfDescendantMaxLinkDistance_ptr(
        std::shared_ptr<tsl::sparse_map<int, IndexPairDistanceStruct>> umaxMaximumOfDescendantMaxLinkDistance) {
        _umaxMaximumOfDescendantMaxLinkDistance = umaxMaximumOfDescendantMaxLinkDistance;
    }

    // VMAX
    inline void set_vmaxMaxCoarseLinkDistance_ptr(std::shared_ptr<tsl::sparse_map<unsigned int, IndexPairDistanceStruct>> umaxMaxCoarseLinkDistance) {
        _vmaxMaxCoarseLinkDistance = umaxMaxCoarseLinkDistance;
    }
    inline void set_vmaxDescendantMaxLinkDistance_ptr(
        std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>>> vmaxDescendantMaxCoarseLinkDistance) {
        _vmaxDescendantMaxCoarseLinkDistance = vmaxDescendantMaxCoarseLinkDistance;
    }
    inline void set_vmaxMaximumOfDescendantMaxLinkDistance_ptr(
        std::shared_ptr<tsl::sparse_map<int, IndexPairDistanceStruct>> vmaxMaximumOfDescendantMaxCoarseLinkDistance) {
        _vmaxMaximumOfDescendantMaxCoarseLinkDistance = vmaxMaximumOfDescendantMaxCoarseLinkDistance;
    }
};

#endif  // PivotLayer_hpp