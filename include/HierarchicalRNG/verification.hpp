// Cole Foster
// 2022-06-03
#include <chrono>
#include <fstream>
#include <vector>

#include "search-sparse-matrix.hpp"
#include "sparse-matrix.hpp"
#include "tsl/sparse_set.h"

namespace verification {

void printSet(tsl::sparse_set<unsigned int> const& set) {
    printf("{,");
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", *it1);
    }
    printf("}\n");
}



void bruteForceRNG(float* const dataPointer, unsigned int const datasetSize, unsigned int const dimension,
                   std::vector<tsl::sparse_set<unsigned int>>& RNG_neighbors, unsigned long long int& distanceComputations, float& time) {
    std::chrono::high_resolution_clock::time_point tStart = std::chrono::high_resolution_clock::now();
    SparseMatrix sparseMatrix(dataPointer, datasetSize, dimension);

    RNG_neighbors.clear();
    RNG_neighbors.resize(datasetSize);

    int statsBatchSize = 100;
    printf("%-10s%-10s%-10s\n", "Q", " ave.Dist", " ave.Time");
    float averageDistance = 0;
    std::chrono::duration<double> averageTime;
    unsigned long long int dStart_query = sparseMatrix.get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_query = std::chrono::high_resolution_clock::now();

    unsigned int index1 = 0, index2 = 0, index3 = 0;
    for (index1 = 0; index1 < datasetSize; index1++) {
        for (index2 = index1 + 1; index2 < datasetSize; index2++) {
            float const distance12 = sparseMatrix._computeDistance(index1, index2);
            bool flag_isNeighbor = true;

            for (index3 = 0; index3 < datasetSize; index3++) {
                if (index3 == index1 || index3 == index2) continue;

                float const distance13 = sparseMatrix._computeDistance(index1, index3);
                if (distance13 < distance12) {
                    float const distance23 = sparseMatrix._computeDistance(index2, index3);
                    if (distance23 < distance12) {
                        flag_isNeighbor = false;
                        break;
                    }
                }
            }

            if (flag_isNeighbor) {
                RNG_neighbors[index1].insert(index2);
                RNG_neighbors[index2].insert(index1);
            }
        }

        // print stats
        if ((index1 + 1) % statsBatchSize == 0) {
            std::chrono::high_resolution_clock::time_point tEnd_query = std::chrono::high_resolution_clock::now();
            unsigned long long int dEnd_query = sparseMatrix.get_distanceComputationCount();
            averageTime = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_query - tStart_query) / statsBatchSize;
            averageDistance = ((float)(dEnd_query - dStart_query)) / statsBatchSize;
            printf("%-10u%-10.2f%-10.3f\n", index1 + 1, averageDistance, averageTime.count() * 1000);
            dStart_query = sparseMatrix.get_distanceComputationCount();
            tStart_query = std::chrono::high_resolution_clock::now();
        }
    }

    std::chrono::high_resolution_clock::time_point tEnd = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();
    distanceComputations = sparseMatrix.get_distanceComputationCount();

    return;
}

/**
 * @brief Get the RNG Neighbors of Queries by Naive approach, O(N^2)
 *
 * @param dimension
 * @param dataPointer
 * @param datasetSize
 * @param testPointer
 * @param RNG_neighbors
 * @param distanceComputations
 * @param time
 */
void bruteForceSearch(unsigned int const dimension, float* const dataPointer, unsigned int const datasetSize, float* const testPointer,
                      unsigned int const testsetSize, std::vector<tsl::sparse_set<unsigned int>>& RNG_neighbors, float& distanceComputations,
                      float& time) {
    std::chrono::high_resolution_clock::time_point tStart = std::chrono::high_resolution_clock::now();
    std::shared_ptr<SparseMatrix> sparseMatrix1 = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension);
    SearchSparseMatrix sparseMatrix(testPointer, testsetSize, dimension, sparseMatrix1);

    RNG_neighbors.clear();
    RNG_neighbors.resize(testsetSize);

    unsigned int queryIndex = 0;
    unsigned int index1 = 0, index2 = 0, index3 = 0;
    for (queryIndex = 0; queryIndex < testsetSize; queryIndex++) {  // for each search query
        unsigned int index1 = datasetSize + queryIndex;
        for (index2 = 0; index2 < datasetSize; index2++) {  // consider all dataset entries as neighbors
            float const distance12 = sparseMatrix._computeDistance(index1, index2);
            bool flag_isNeighbor = true;

            for (index3 = 0; index3 < datasetSize; index3++) {
                if (index3 == index1 || index3 == index2) continue;

                float const distance13 = sparseMatrix._computeDistance(index1, index3);
                if (distance13 < distance12) {
                    float const distance23 = sparseMatrix._computeDistance(index2, index3);
                    if (distance23 < distance12) {
                        flag_isNeighbor = false;
                        break;
                    }
                }
            }

            if (flag_isNeighbor) {
                RNG_neighbors[queryIndex].insert(index2);
            }
        }
    }

    std::chrono::high_resolution_clock::time_point tEnd = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() * 1000 / testsetSize;
    distanceComputations = (float)(((float)sparseMatrix.get_distanceComputationCount()) / ((float)testsetSize));

    return;
}

// compare G2 to G1 as the ground truth graph. return true if exactly the same.
bool graphComparison(std::vector<tsl::sparse_set<unsigned int>>& G_GT, std::vector<tsl::sparse_set<unsigned int>>& G, unsigned int& numberOfLinks,
                     unsigned int& extraLinks, unsigned int& missingLinks) {
    bool success = true;
    numberOfLinks = 0;
    extraLinks = 0;
    missingLinks = 0;

    // ensure correct number of indices
    unsigned int numberOfIndices = (unsigned int)G_GT.size();
    if (G.size() != numberOfIndices) {
        success = false;
        printf("Incorrect number of indices given...\n");
        return success;
    }

    // check all neighbors for correctness
    for (unsigned int index1 = 0; index1 < numberOfIndices; index1++) {
        tsl::sparse_set<unsigned int> const& neighbors1_GT = G_GT[index1];
        tsl::sparse_set<unsigned int> const& neighbors1 = G[index1];
        numberOfLinks += neighbors1_GT.size();

        // see if correct
        if (neighbors1_GT != neighbors1) {
            success = false;

            // find missing links
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = neighbors1_GT.begin(); it2 != neighbors1_GT.end(); it2++) {
                unsigned int index2 = (*it2);
                if (neighbors1.find(index2) == neighbors1.end()) {
                    missingLinks += 1;
                }
            }

            // find extra links
            for (it2 = neighbors1.begin(); it2 != neighbors1.end(); it2++) {
                unsigned int index2 = (*it2);
                if (neighbors1_GT.find(index2) == neighbors1_GT.end()) {
                    extraLinks += 1;
                }
            }
        }
    }
    numberOfLinks /= 2;
    missingLinks /= 2;
    extraLinks /= 2;

    return success;
}

// compare G2 to G1 as the ground truth graph. return true if exactly the same.
bool neighborsComparison(std::vector<tsl::sparse_set<unsigned int>>& G_GT, std::vector<tsl::sparse_set<unsigned int>>& G, unsigned int& numberOfLinks,
                         unsigned int& extraLinks, unsigned int& missingLinks) {
    bool success = true;
    numberOfLinks = 0;
    extraLinks = 0;
    missingLinks = 0;

    // ensure correct number of indices
    unsigned int numberOfIndices = (unsigned int)G_GT.size();
    if (G.size() != numberOfIndices) {
        success = false;
        printf("Incorrect number of indices given...\n");
        return success;
    }

    // check all neighbors for correctness
    for (unsigned int index1 = 0; index1 < numberOfIndices; index1++) {
        tsl::sparse_set<unsigned int> const& neighbors1_GT = G_GT[index1];
        tsl::sparse_set<unsigned int> const& neighbors1 = G[index1];
        numberOfLinks += neighbors1_GT.size();

        // printSet(neighbors1_GT);
        // printSet(neighbors1);

        // see if correct
        if (neighbors1_GT != neighbors1) {
            success = false;

            // find missing links
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = neighbors1_GT.begin(); it2 != neighbors1_GT.end(); it2++) {
                unsigned int index2 = (*it2);
                if (neighbors1.find(index2) == neighbors1.end()) {
                    missingLinks += 1;
                }
            }

            // find extra links
            for (it2 = neighbors1.begin(); it2 != neighbors1.end(); it2++) {
                unsigned int index2 = (*it2);
                if (neighbors1_GT.find(index2) == neighbors1_GT.end()) {
                    extraLinks += 1;
                }
            }
        }
    }

    return success;
}

}  // namespace verification