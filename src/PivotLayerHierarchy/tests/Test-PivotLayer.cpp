#include "Test-PivotLayer.hpp"

/*
    void test_constructors();
    void test_addPivotToLayer();
    void test_initializePivotsInLayer();
    void test_sharedPointers();
    void test_pivotLayerNeighbors();
    void test_sparseMatrix();
    void test_umaxMaxLinkDistance();
    void test_ancestors();
    void test_descendants();
    void test_bruteForceConstruction();
    void test_copyConstructor();
*/

void Test_PivotLayer::test_constructors() {
    // default constructor
    PivotLayer layer1;
    BOOST_CHECK(layer1.get_pivotLayerID() == 0);
    BOOST_CHECK(layer1.get_numberOfPivotLayers() == 1);
    BOOST_CHECK(layer1.get_luneType() == "Berk");

    // normal constructor
    std::shared_ptr<SparseMatrix> matrix;
    PivotLayer layer2(2, 7, matrix, "Cole", true);
    BOOST_CHECK(layer2.get_pivotLayerID() == 2);
    BOOST_CHECK(layer2.get_numberOfPivotLayers() == 7);
    BOOST_CHECK(layer2.get_luneType() == "Cole");
}

void Test_PivotLayer::test_addPivotToLayer() {
    // intialize layer
    std::shared_ptr<SparseMatrix> matrix;
    PivotLayer layer(2, 5, matrix, "Berk", false);

    // add random things to layer
    tsl::sparse_set<unsigned int> indices;
    tsl::sparse_map<unsigned int, float> radii;
    for (unsigned int i = 0; i < 200; i++) {
        float rad = 3.1415926 * (float)i;
        layer.addPivotToLayer(i, rad);
        indices.insert(i);
        radii.emplace(i, rad);
    }

    // check we get the correct results!
    BOOST_CHECK(indices == *layer.get_pivotIndices_ptr());
    BOOST_CHECK(radii == *layer.get_pivotRadii_ptr());
}

void Test_PivotLayer::test_initializePivotsInLayer() {
    // initialize top, middle, bottom layers
    int num_layers = 3;
    std::shared_ptr<SparseMatrix> matrix;
    PivotLayer layer1(0, num_layers, matrix, "Berk", false);
    PivotLayer layer2(1, num_layers, matrix, "Berk", false);
    PivotLayer layer3(2, num_layers, matrix, "Berk", false);

    // initialize some structures
    unsigned int num_pivots = 100;
    for (unsigned int i = 0; i < num_pivots; i++) {
        layer1.initializePivotInLayer(i);
        layer2.initializePivotInLayer(i);
        layer3.initializePivotInLayer(i);
    }

    // top layer
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& neighbors1 = *layer1.get_pivotLayerNeighbors_ptr();
    BOOST_CHECK(neighbors1.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& coarse1 = *layer1.get_coarsePivotLayerNeighbors_ptr();
    BOOST_CHECK(coarse1.size() == 0);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> const& ancestors1 = *layer1.get_ancestorPivotIndices_ptr();
    BOOST_CHECK(ancestors1.size() == 0);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> const& descendants1 = *layer1.get_descendantPivotIndices_ptr();
    BOOST_CHECK(descendants1.size() == num_pivots);
    BOOST_CHECK(descendants1.at(3).size() == 2);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>> const& dmax1 = *layer1.get_dmaxMaxChildDistance_ptr();
    BOOST_CHECK(dmax1.size() == num_pivots);
    BOOST_CHECK(dmax1.at(3).size() == 2);
    tsl::sparse_map<unsigned int, IndexPairDistanceStruct> const& umax1 = *layer1.get_umaxMaxLinkDistance_ptr();
    BOOST_CHECK(umax1.size() == 0);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>> const& umaxDescendants1 = *layer1.get_umaxDescendantMaxLinkDistance_ptr();
    BOOST_CHECK(umaxDescendants1.size() == num_pivots);
    BOOST_CHECK(umaxDescendants1.at(3).size() == 2);

    // middle layer
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& neighbors2 = *layer2.get_pivotLayerNeighbors_ptr();
    BOOST_CHECK(neighbors2.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& coarse2 = *layer2.get_coarsePivotLayerNeighbors_ptr();
    BOOST_CHECK(coarse2.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> const& ancestors2 = *layer2.get_ancestorPivotIndices_ptr();
    BOOST_CHECK(ancestors2.size() == num_pivots);
    BOOST_CHECK(ancestors2.at(3).size() == 1);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> const& descendants2 = *layer2.get_descendantPivotIndices_ptr();
    BOOST_CHECK(descendants2.size() == num_pivots);
    BOOST_CHECK(descendants2.at(3).size() == 1);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>> const& dmax2 = *layer2.get_dmaxMaxChildDistance_ptr();
    BOOST_CHECK(dmax2.size() == num_pivots);
    BOOST_CHECK(dmax2.at(3).size() == 1);
    tsl::sparse_map<unsigned int, IndexPairDistanceStruct> const& umax2 = *layer2.get_umaxMaxLinkDistance_ptr();
    BOOST_CHECK(umax2.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>> const& umaxDescendants2 = *layer2.get_umaxDescendantMaxLinkDistance_ptr();
    BOOST_CHECK(umaxDescendants2.size() == num_pivots);
    BOOST_CHECK(umaxDescendants2.at(3).size() == 1);

    // bottom layer
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& neighbors3 = *layer3.get_pivotLayerNeighbors_ptr();
    BOOST_CHECK(neighbors3.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& coarse3 = *layer3.get_coarsePivotLayerNeighbors_ptr();
    BOOST_CHECK(coarse3.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> const& ancestors3 = *layer3.get_ancestorPivotIndices_ptr();
    BOOST_CHECK(ancestors3.size() == num_pivots);
    BOOST_CHECK(ancestors3.at(3).size() == 2);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, tsl::sparse_set<unsigned int>>> const& descendants3 = *layer3.get_descendantPivotIndices_ptr();
    BOOST_CHECK(descendants3.size() == 0);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, float>> const& dmax3 = *layer3.get_dmaxMaxChildDistance_ptr();
    BOOST_CHECK(dmax3.size() == 0);
    tsl::sparse_map<unsigned int, IndexPairDistanceStruct> const& umax3 = *layer3.get_umaxMaxLinkDistance_ptr();
    BOOST_CHECK(umax3.size() == num_pivots);
    tsl::sparse_map<unsigned int, tsl::sparse_map<int, IndexPairDistanceStruct>> const& umaxDescendants3 = *layer3.get_umaxDescendantMaxLinkDistance_ptr();
    BOOST_CHECK(umaxDescendants3.size() == 0);
}

void Test_PivotLayer::test_sharedPointers() {
    PivotLayer pivotLayer;

    // create indices set
    std::shared_ptr<tsl::sparse_set<unsigned int>> indices1 = std::make_shared<tsl::sparse_set<unsigned int>>();
    for (unsigned int i = 0; i < 100; i++) {
        (*indices1).insert(i);
    }

    // set for the pivot layer
    pivotLayer.set_pivotIndices_ptr(indices1);
    (*indices1).insert({101, 102, 103});  // adding pivots after setting

    // these should be the same, otherwise we are just making a copy!
    tsl::sparse_set<unsigned int> const& indices2 = *pivotLayer.get_pivotIndices_ptr();
    BOOST_CHECK((*indices1) == indices2);
}

void Test_PivotLayer::test_sparseMatrix() {
    bool success = true;

    float data[10] = {6, 8, 3, 4, 3, 6, 4, 3, 4, 3};
    float* dataPointer = data;
    std::shared_ptr<SparseMatrix> SM = std::make_shared<SparseMatrix>(dataPointer, 5, 2);
    PivotLayer layer(2, 4, SM, "Berk", false);

    // ensure distance computation count starts at zero
    BOOST_CHECK(layer.get_distanceComputationCount() == 0);

    // check correct computation
    float const distance = layer.getDistance(0, 1);
    BOOST_CHECK(distance == (float)5.0f);
    BOOST_CHECK(layer.get_distanceComputationCount() == 1);

    // get distance again
    float const distance1 = layer.getDistance(0, 1);
    BOOST_CHECK(distance1 == (float)5.0f);
    BOOST_CHECK(layer.get_distanceComputationCount() == 1);

    // get distance other side
    float const distance2 = layer.getDistance(1, 0);
    BOOST_CHECK(distance2 == (float)5.0f);
    BOOST_CHECK(layer.get_distanceComputationCount() == 1);

    // get new distance, but same index
    // shouldn't compute distance, just return 0
    float distance3 = layer.getDistance(2, 2);
    BOOST_CHECK(distance3 == 0.0f);
    BOOST_CHECK(layer.get_distanceComputationCount() == 1);
}

void Test_PivotLayer::test_copyConstructor() {
    // intialize layer
    std::shared_ptr<SparseMatrix> matrix;
    PivotLayer layer1(2, 5, matrix, "Berk", false);

    // add random things to layer
    for (unsigned int i = 0; i < 10; i++) {
        layer1.addPivotToLayer(i, 3.1415926 * (float)i);
        layer1.initializePivotInLayer(i);
    }

    // copy constructor
    PivotLayer layer2(layer1);

    // check equality of all structures, should all be equal since sharing pointers!
    BOOST_CHECK(layer1.get_pivotLayerID() == layer2.get_pivotLayerID());
    BOOST_CHECK(layer1.get_numberOfPivotLayers() == layer2.get_numberOfPivotLayers());
    BOOST_CHECK(layer1.get_luneType() == layer2.get_luneType());
    BOOST_CHECK(layer1.get_pivotIndices_ptr() == layer2.get_pivotIndices_ptr());
    BOOST_CHECK(layer1.get_pivotRadii_ptr() == layer2.get_pivotRadii_ptr());
    BOOST_CHECK(layer1.get_pivotLayerNeighbors_ptr() == layer2.get_pivotLayerNeighbors_ptr());
    BOOST_CHECK(layer1.get_coarsePivotLayerNeighbors_ptr() == layer2.get_coarsePivotLayerNeighbors_ptr());
    BOOST_CHECK(layer1.get_ancestorPivotIndices_ptr() == layer2.get_ancestorPivotIndices_ptr());
    BOOST_CHECK(layer1.get_descendantPivotIndices_ptr() == layer2.get_descendantPivotIndices_ptr());
    BOOST_CHECK(layer1.get_dmaxMaxChildDistance_ptr() == layer2.get_dmaxMaxChildDistance_ptr());
    BOOST_CHECK(layer1.get_dmaxMaximumOfMaxChildDistance_ptr() == layer2.get_dmaxMaximumOfMaxChildDistance_ptr());
    BOOST_CHECK(layer1.get_umaxMaxLinkDistance_ptr() == layer2.get_umaxMaxLinkDistance_ptr());
    BOOST_CHECK(layer1.get_umaxDescendantMaxLinkDistance_ptr() == layer2.get_umaxDescendantMaxLinkDistance_ptr());
    BOOST_CHECK(layer1.get_umaxMaximumOfDescendantMaxLinkDistance_ptr() == layer2.get_umaxMaximumOfDescendantMaxLinkDistance_ptr());

    // add some stuff to layer1
    for (unsigned int i = 10; i < 30; i++) {
        layer1.addPivotToLayer(i, 3.1415926 * (float)i);
        layer1.initializePivotInLayer(i);
    }

    // check equality of all structures, should all be equal since sharing pointers!
    BOOST_CHECK(layer1.get_pivotLayerID() == layer2.get_pivotLayerID());
    BOOST_CHECK(layer1.get_numberOfPivotLayers() == layer2.get_numberOfPivotLayers());
    BOOST_CHECK(layer1.get_luneType() == layer2.get_luneType());
    BOOST_CHECK(layer1.get_pivotIndices_ptr() == layer2.get_pivotIndices_ptr());
    BOOST_CHECK(layer1.get_pivotRadii_ptr() == layer2.get_pivotRadii_ptr());
    BOOST_CHECK(layer1.get_pivotLayerNeighbors_ptr() == layer2.get_pivotLayerNeighbors_ptr());
    BOOST_CHECK(layer1.get_coarsePivotLayerNeighbors_ptr() == layer2.get_coarsePivotLayerNeighbors_ptr());
    BOOST_CHECK(layer1.get_ancestorPivotIndices_ptr() == layer2.get_ancestorPivotIndices_ptr());
    BOOST_CHECK(layer1.get_descendantPivotIndices_ptr() == layer2.get_descendantPivotIndices_ptr());
    BOOST_CHECK(layer1.get_dmaxMaxChildDistance_ptr() == layer2.get_dmaxMaxChildDistance_ptr());
    BOOST_CHECK(layer1.get_dmaxMaximumOfMaxChildDistance_ptr() == layer2.get_dmaxMaximumOfMaxChildDistance_ptr());
    BOOST_CHECK(layer1.get_umaxMaxLinkDistance_ptr() == layer2.get_umaxMaxLinkDistance_ptr());
    BOOST_CHECK(layer1.get_umaxDescendantMaxLinkDistance_ptr() == layer2.get_umaxDescendantMaxLinkDistance_ptr());
    BOOST_CHECK(layer1.get_umaxMaximumOfDescendantMaxLinkDistance_ptr() == layer2.get_umaxMaximumOfDescendantMaxLinkDistance_ptr());
}

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

void printMapSet(tsl::sparse_map<unsigned, tsl::sparse_set<unsigned int>> mapSet) {
    std::cout << "{";
    tsl::sparse_map<unsigned, tsl::sparse_set<unsigned int>>::const_iterator it1;
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

void Test_PivotLayer::test_lune() {
    std::shared_ptr<SparseMatrix> SM;
    PivotLayer layer1(0, 2, SM, "Berk", false);
    PivotLayer layer2(0, 2, SM, "Cole", false);

    float distance = 2.342;
    float radius1 = 0.234;
    float radius2 = 0.423;

    BOOST_CHECK(layer1._luneRadius(distance, radius1, radius2) == (distance - 2 * radius1 - radius2));
    BOOST_CHECK(layer1._luneRadius(distance, radius2, radius1) == (distance - 2 * radius2 - radius1));
    BOOST_CHECK(layer2._luneRadius(distance, radius1, radius2) == (distance - radius2));
    BOOST_CHECK(layer2._luneRadius(distance, radius2, radius1) == (distance - radius1));
}

void Test_PivotLayer::test_bruteForceConstruction() {
    float* dataPointer;
    unsigned int D = 7;
    unsigned int N = 100;
    uniformDataset(dataPointer, D, N, 3);

    std::shared_ptr<SparseMatrix> SM = std::make_shared<SparseMatrix>(dataPointer, N, D);
    PivotLayer layer(0, 3, SM, "Berk", false);

    // compute brute force berk graph
    float radius = 0.2343;
    for (unsigned int i = 0; i < N; i++) {
        layer.addPivotToLayer(i, radius);
        layer.initializePivotInLayer(i);
    }
    layer.computeBruteForcePivotLayerGraph();

    // now, do it by brute force here
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> BG;
    for (unsigned int i = 0; i < N; i++) {
        BG.emplace(i, tsl::sparse_set<unsigned int>{});
    }

    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
            if (j == i) continue;
            float distance12 = SM->getDistance(i, j);
            bool isNeighbors = true;
            for (unsigned int k = 0; k < N; k++) {
                if (k == i || k == j) continue;
                float distance13 = SM->getDistance(i, k);
                if (distance13 < distance12 - 3 * radius) {
                    float distance23 = SM->getDistance(j, k);
                    if (distance23 < distance12 - 3 * radius) {
                        isNeighbors = false;
                    }
                }
            }

            if (isNeighbors) {
                BG[i].insert(j);
                BG[j].insert(i);
            }
        }
    }

    BOOST_CHECK(BG == *layer.get_pivotLayerNeighbors_ptr());

    delete[] dataPointer;
}