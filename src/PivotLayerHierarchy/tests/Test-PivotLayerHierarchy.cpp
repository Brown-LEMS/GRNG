#include "Test-PivotLayerHierarchy.hpp"

// generate dataset of [N,D] floating point numbers
// uniformly distributed from [-1,1]
void uniformDataset(float*& dataPointer, unsigned int dimension, unsigned int N, unsigned int seed) {
    srand(seed);  // for same initializations of random variables

    unsigned int numElements = (unsigned int)dimension * N;
    // delete[] dataPointer;
    dataPointer = new float[N * dimension];

    for (int i = 0; i < numElements; i++) {
        dataPointer[i] = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * 2 - 1;
    }
}

void Test_PivotLayerHierarchy::test_constructors() {
    // test default constructor
    PivotLayerHierarchy pivotLayers1;
    BOOST_CHECK(pivotLayers1.get_pivotRadiusVector() == std::vector<float>{0.0f});
    BOOST_CHECK(pivotLayers1.get_numberOfPivotLayers() == 1);
    BOOST_CHECK(pivotLayers1.get_luneType() == "Berk");

    // test other constructor
    std::vector<float> radiiVec = {0.234, 0.3432, 0};
    unsigned int datasetSize = 100;
    unsigned int dimension = 3;
    float* dataPointer;
    uniformDataset(dataPointer, dimension, datasetSize, 3);
    PivotLayerHierarchy pivotLayers3(radiiVec, "Cole", dataPointer, datasetSize, dimension);
    BOOST_CHECK(pivotLayers3.get_pivotRadiusVector() == radiiVec);
    BOOST_CHECK(pivotLayers3.get_numberOfPivotLayers() == radiiVec.size());
    BOOST_CHECK(pivotLayers3.get_luneType() == "Cole");
    BOOST_CHECK(pivotLayers3.get_pivotLayers().size() == 3);

    delete[] dataPointer;
}

void Test_PivotLayerHierarchy::test_copyConstructor() {
    unsigned int datasetSize = 100;
    unsigned int dimension = 3;
    float* dataPointer;
    uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::vector<float> radiiVec = {0.312, 0.123, 0};
    PivotLayerHierarchy pivotLayers(radiiVec, "Cole", dataPointer, datasetSize, dimension);
    float d1 = pivotLayers.getDistance(0, 1);
    float d2 = pivotLayers.getDistance(1, 4);
    float d3 = pivotLayers.getDistance(3, 45);

    // copy constructor
    PivotLayerHierarchy pivotLayers1(pivotLayers);
    BOOST_CHECK(pivotLayers1.get_pivotRadiusVector() == radiiVec);
    BOOST_CHECK(pivotLayers1.get_numberOfPivotLayers() == radiiVec.size());
    BOOST_CHECK(pivotLayers1.get_luneType() == "Cole");
    BOOST_CHECK(pivotLayers1.get_pivotLayers().size() == 3);

    // check same SparseMatrix and pivot layers
    BOOST_CHECK(pivotLayers1.get_distanceComputationCount() == 3);
    BOOST_CHECK(pivotLayers1.get_distanceComputationCount() == pivotLayers.get_distanceComputationCount());

    BOOST_CHECK(pivotLayers1.get_pivotLayers().size() == pivotLayers.get_pivotLayers().size());
}

void Test_PivotLayerHierarchy::test_sparseMatrix() {
    unsigned int datasetSize = 5;
    unsigned int dimension = 2;
    float data[10] = {6, 8, 3, 4, 3, 6, 4, 3, 4, 3};
    float* dataPointer = data;
    PivotLayerHierarchy pivotLayers({0}, "Cole", dataPointer, datasetSize, dimension);

    BOOST_CHECK(pivotLayers.get_distanceComputationCount() == 0);

    SparseMatrix& SM = *pivotLayers.get_sparseMatrix();
    BOOST_CHECK(SM.getDistance(0, 1) == 5.0f);
    BOOST_CHECK(SM.get_distanceComputationCount() == 1);
    BOOST_CHECK(pivotLayers.get_distanceComputationCount() == 1);

    // get distance again
    float const distance1 = pivotLayers.getDistance(0, 1);
    BOOST_CHECK(distance1 == (float)5.0f);
    BOOST_CHECK(pivotLayers.get_distanceComputationCount() == 1);

    // get distance other side
    float const distance2 = pivotLayers.getDistance(1, 0);
    BOOST_CHECK(distance2 == (float)5.0f);
    BOOST_CHECK(pivotLayers.get_distanceComputationCount() == 1);

    // get new distance, but same index
    // shouldn't compute distance, just return 0
    float distance3 = pivotLayers.getDistance(2, 2);
    BOOST_CHECK(distance3 == 0.0f);
    BOOST_CHECK(pivotLayers.get_distanceComputationCount() == 1);
}

void print(tsl::sparse_set<unsigned int> const& set) {
    std::cout << "{";
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        std::cout << (*it1) << ", ";
    }
    std::cout << "} " << std::endl;
};

void print(tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& mapSet) {
    std::cout << "{";
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = mapSet.begin(); it1 != mapSet.end(); it1++) {
        unsigned int key = (*it1).first;
        std::cout << key << ": {";
        tsl::sparse_set<unsigned int> set = (*it1).second;
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = set.begin(); it1 != set.end(); it1++) {
            std::cout << (*it1) << ", ";
        }
        std::cout << "} ";
    }
    std::cout << "}" << std::endl;
};

void Test_PivotLayerHierarchy::test_ancestryDotCom() {
    unsigned int datasetSize = 100;
    unsigned int dimension = 3;
    float* dataPointer;
    uniformDataset(dataPointer, dimension, datasetSize, 3);
    PivotLayerHierarchy pivotLayers({0.312, 0.12, 0}, "Cole", dataPointer, datasetSize, dimension);

    // add a bunch to top layer and bottom layer
    for (unsigned int i = 2; i < datasetSize; i += 3) {
        pivotLayers.addPivotToLayer(0, i - 2, 0.312);
        pivotLayers.addPivotToLayer(1, i - 1, 0.12);
        pivotLayers.addPivotToLayer(2, i, 0);
    }

    // add some random ancestor relationships
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> ancestors1;
    ancestors1.emplace(0, tsl::sparse_set<unsigned int>{0, 3, 6, 9, 12, 15, 18, 21});
    ancestors1.emplace(1, tsl::sparse_set<unsigned int>{1, 4, 7, 10});
    pivotLayers.updateAncestors(2, 5, ancestors1);
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> ancestors2;
    ancestors2.emplace(0, tsl::sparse_set<unsigned int>{12, 15, 18, 21});
    ancestors2.emplace(1, tsl::sparse_set<unsigned int>{1, 10});
    pivotLayers.updateAncestors(2, 11, ancestors2);
    tsl::sparse_map<int, tsl::sparse_set<unsigned int>> ancestors3;
    ancestors1.emplace(0, tsl::sparse_set<unsigned int>{0, 3, 6, 9, 12, 15, 18, 21, 33, 96});
    ancestors1.emplace(1, tsl::sparse_set<unsigned int>{1, 4, 7, 10, 13, 14});
    pivotLayers.updateAncestors(2, 8, ancestors3);

    // now, check that ancestors and descendants are correct. just the bottom layer
    int const layerIndex = 2;
    tsl::sparse_set<unsigned int> const& pivotIndices = *pivotLayers[layerIndex].get_pivotIndices_ptr();
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>> const& ancestorIndicesMap =
            pivotLayers[layerIndex].get_ancestorPivotIndicesMap(pivotIndex);

        // std::cout << pivotIndex << " parents: " << std::endl;
        // print(ancestorIndicesMap);

        // for each above layer
        tsl::sparse_map<int, tsl::sparse_set<unsigned int>>::const_iterator it2;
        for (it2 = ancestorIndicesMap.begin(); it2 != ancestorIndicesMap.end(); it2++) {
            int const ancestorLayerIndex = (*it2).first;
            tsl::sparse_set<unsigned int> const& ancestorIndices = (*it2).second;

            // for each ancestor
            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = ancestorIndices.begin(); it3 != ancestorIndices.end(); it3++) {
                unsigned int const ancestorIndex = (*it3);
                tsl::sparse_set<unsigned int> const& descendants =
                    pivotLayers[ancestorLayerIndex].get_descendantPivotIndices(ancestorIndex, layerIndex);

                // std::cout << "Layer " << ancestorLayerIndex << "-" << ancestorIndex << " children: " << std::endl;
                // print(descendants);

                // check that parent/child relationship matches
                BOOST_REQUIRE(descendants.find(pivotIndex) != descendants.end());  // in
            }
        }
    }

    delete[] dataPointer;
    return;
}

void Test_PivotLayerHierarchy::test_umax() {
    unsigned int datasetSize = 100;
    unsigned int dimension = 3;
    float* dataPointer;
    uniformDataset(dataPointer, dimension, datasetSize, 3);
    PivotLayerHierarchy pivotLayers({0.312, 0.12, 0}, "Berk", dataPointer, datasetSize, dimension);
    int const layerIndex = 2;

    // add a bunch to top layer and bottom layer
    for (unsigned int i = 0; i < datasetSize; i++) {
        pivotLayers.addPivotToLayer(layerIndex, i, 0.011);
    }

    // add some links in the layer
    pivotLayers.addPivotLayerNeighbors(layerIndex, 0, tsl::sparse_set<unsigned int>{3, 6, 7, 4, 10, 23, 44, 5, 43, 2});
    pivotLayers.addPivotLayerNeighbors(layerIndex, 10, tsl::sparse_set<unsigned int>{23, 44, 5, 43, 2});
    pivotLayers.addPivotLayerNeighbors(layerIndex, 66, tsl::sparse_set<unsigned int>{3, 6, 7, 4, 10, 45, 99, 23, 44, 5, 43, 2});
    pivotLayers.addPivotLayerNeighbors(layerIndex, 23, tsl::sparse_set<unsigned int>{3, 6, 7, 4, 10, 52, 5, 43, 2});

    // for each pivot in the layer
    tsl::sparse_set<unsigned int> const& pivotIndices = *pivotLayers[layerIndex].get_pivotIndices_ptr();
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        float const pivotRadius = pivotLayers[layerIndex].get_pivotRadius(pivotIndex);
        IndexPairDistanceStruct const& umax = pivotLayers[layerIndex].get_umaxMaxLinkDistance(pivotIndex);

        // find max link distance
        tsl::sparse_set<unsigned int> const& neighbors = pivotLayers[layerIndex].get_pivotLayerNeighbors(pivotIndex);
        IndexPairDistanceStruct umax_pivot;

        // find link distance for each neighbor
        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = neighbors.begin(); it2 != neighbors.end(); it2++) {
            int const neighborIndex = (*it2);
            float const neighborRadius = pivotLayers[layerIndex].get_pivotRadius(neighborIndex);
            float const link_distance = pivotLayers.getDistance(pivotIndex, neighborIndex) - (2 * pivotRadius + neighborRadius);
            IndexPairDistanceStruct umax_temp(pivotIndex, neighborIndex, link_distance);

            if (umax_pivot.isInvalid()) {
                umax_pivot = umax_temp;
            } else if (link_distance > umax_pivot.distance) {
                umax_pivot = umax_temp;
            }
        }

        if (umax != umax_pivot) {
            std::cout << "real: [" << umax.index1 << "," << umax.index2 << "]" << std::endl;
            std::cout << "thought: [" << umax_pivot.index1 << "," << umax_pivot.index2 << "]" << std::endl;
            BOOST_FAIL("incorrect umax");
        }
    }

    delete[] dataPointer;
    return;
}

void Test_PivotLayerHierarchy::test_umaxDescendants(){};

void Test_PivotLayerHierarchy::test_umaxMaximumOfDescendants(){};
