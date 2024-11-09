#include "Test_OnlineSearch.hpp"

#include "datasets.hpp"
#include "incremental-build.hpp"
#include "offline-build.hpp"

void Test_OnlineSearch::test_neighbors() {
    std::cout << "Unit Test: Neighbors" << std::endl;

    // create dataset
    unsigned int const datasetSize = 1600;  // 3200;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);

    // 3L hierarchy, 2D, N=void test_buildAndDataStructures();
    // std::vector<float> pivotRadiusVector = std::vector<float>{0.21872, 0};
    // std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0};
    std::vector<float> pivotRadiusVector = std::vector<float>{0.21872, 0.08995, 0.02123, 0};
    std::vector<float> effectiveRadiusVector = Test_OnlineSearch::getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension,true);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.setNumThreads(1);
    IB.build(false);

    // create testset
    unsigned int const testsetSize = 100;  // 3200;
    float* testPointer;
    Datasets::uniformDataset(testPointer, dimension, testsetSize, 5);

    // perform online search
    OnlineSearch OS(pivotLayers);
    OS.setNumThreads(8);
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> searchNeighbors_hRNG{};
    OS.search(testPointer, testsetSize, dimension, true, searchNeighbors_hRNG);

    // perform brute force search
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> searchNeighbors_GT{};
    getRNGNeighbors(dataPointer, datasetSize, dimension, testPointer, testsetSize, searchNeighbors_GT);

    // test for correct RNG neighbors
    bool flag_success = true;
    for (int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        if (searchNeighbors_GT[queryIndex] != searchNeighbors_hRNG[queryIndex]) {
            printf("Incorrect Search RNG for query %u: \n",queryIndex);
            printf("GT: "); printSet(searchNeighbors_GT[queryIndex]);
            printf("hRNG: "); printSet(searchNeighbors_hRNG[queryIndex]);
            flag_success = false;
            break;
        }
    }

    if (!flag_success) {
        delete[] dataPointer;
        delete[] testPointer;       
        BOOST_FAIL("Incorrect Search Neighbors...");
    }

    delete[] dataPointer;
    delete[] testPointer;
    return;
}

void Test_OnlineSearch::getRNGNeighbors(float* dataPointer, unsigned int const datasetSize, unsigned int const dimension, float* testPointer,
                     unsigned int const testsetSize, tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& RNG) {
    std::shared_ptr<SparseMatrix> sparseMatrix_data = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension);
    SearchSparseMatrix sparseMatrix(testPointer, testsetSize, dimension, sparseMatrix_data);

    RNG.clear();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        RNG[queryIndex] = tsl::sparse_set<unsigned int>{};
    }

    // consider each query
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        unsigned int const index1 = datasetSize + queryIndex;

        // consider all indices as a neighbor
        for (int index2 = 0; index2 < datasetSize; index2++) {
            bool flag_isNeighbor = true;
            float const distance12 = sparseMatrix.getDistance(index1, index2);

            // check for interference by other pivots
            for (int index3 = 0; index3 < datasetSize; index3++) {
                if (index3 == index2) continue;

                // check both sides for interference
                float const distance13 = sparseMatrix.getDistance(index1, index3);
                if (distance13 < distance12) {
                    float const distance23 = sparseMatrix.getDistance(index2, index3);
                    if (distance23 < distance12) {
                        flag_isNeighbor = false;
                        break;
                    }
                }
            }

            if (flag_isNeighbor) {
                RNG[queryIndex].insert(index2);
                continue;
            }
        }
    }

    printf("Done with Neighbors test");
    return;
}