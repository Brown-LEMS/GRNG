// Cole Foster
// 2022-06-03

// Implementation of ApproximateRNG by Hacid et. al
// Hacid, Hakim, and Tetsuya Yoshida. 2007. “Incremental Neighborhood Graphs Construction for Multidimensional Databases Indexing.” In Advances in
// Artificial Intelligence, edited by Ziad Kobti and Dan Wu, 405–16. Lecture Notes in Computer Science. Berlin, Heidelberg: Springer.
// https://doi.org/10.1007/978-3-540-72665-4_35.
#include <omp.h>

#include <chrono>
#include <vector>

#include "links-map.hpp"
#include "search-sparse-matrix.hpp"
#include "sparse-matrix.hpp"
#include "tsl/sparse_set.h"

namespace hacid {

// naive nearest neighbor search
void nearest_neighbor(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix& sparseMatrix, int numThreads,
                      unsigned int& nearestNeighbor, float& nearestDistance) {
    nearestDistance = HUGE_VAL;

#pragma omp parallel num_threads(numThreads)
    {
        SparseMatrixLocal localSparseMatrix = sparseMatrix.getLocalCopy();
        unsigned int thread_nn = 0;
        float thread_nd = HUGE_VAL;

#pragma omp for schedule(static)
        for (unsigned int candidateIndex = 0; candidateIndex < datasetSize; candidateIndex++) {
            float const distance = localSparseMatrix._computeDistance(queryIndex, candidateIndex);
            if (distance < thread_nd) {
                thread_nd = distance;
                thread_nn = candidateIndex;
            }
        }

#pragma omp critical(updateNN)
        {
            if (thread_nd < nearestDistance) {
                nearestDistance = thread_nd;
                nearestNeighbor = thread_nn;
            }
            sparseMatrix.updateSparseMatrixWithLocal(localSparseMatrix);
        }
    }

    return;
}

// finding the furthest distance to a neighbor
float furthest_neighbor(unsigned int const index1, SparseMatrix& sparseMatrix, std::vector<tsl::sparse_set<unsigned int>> const& RNG_neighbors) {
    float furthestDistance = 0;

    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = RNG_neighbors[index1].begin(); it1 != RNG_neighbors[index1].end(); it1++) {
        unsigned int const index2 = (*it1);
        float const distance = sparseMatrix._computeDistance(index1, index2);
        if (distance > furthestDistance) {
            furthestDistance = distance;
        }
    }

    return furthestDistance;
}

// naive nearest neighbor search
void collectSRIndices(unsigned int const queryIndex, unsigned int const datasetSize, float const SR, SparseMatrix& sparseMatrix, int numThreads,
                      tsl::sparse_set<unsigned int>& SR_region) {
    SR_region.clear();

#pragma omp parallel num_threads(numThreads)
    {
        tsl::sparse_set<unsigned int> thread_SR{};
        SparseMatrixLocal localSparseMatrix = sparseMatrix.getLocalCopy();

#pragma omp for schedule(static)
        for (unsigned int candidateIndex = 0; candidateIndex < datasetSize; candidateIndex++) {
            float const distance = localSparseMatrix._computeDistance(queryIndex, candidateIndex);

            if (distance <= SR) {
                thread_SR.insert(candidateIndex);
            }
        }

#pragma omp critical(updateSR)
        { 
            SR_region.insert(thread_SR.begin(), thread_SR.end()); 
        sparseMatrix.updateSparseMatrixWithLocal(localSparseMatrix);
        }
    }

    return;
}

// get RNG neighbors of query using indices within SR_region
void getRNGNeighbors(unsigned int const queryIndex, tsl::sparse_set<unsigned int> const& SR_region, SparseMatrix& sparseMatrix,
                     std::vector<tsl::sparse_set<unsigned int>>& RNG_neighbors) {
    RNG_neighbors[queryIndex].clear();

    // consider all within region as interfering
    tsl::sparse_set<unsigned int>::const_iterator it2, it3;
    for (it2 = SR_region.begin(); it2 != SR_region.end(); it2++) {
        unsigned int const index2 = (*it2);
        float const distance12 = sparseMatrix._computeDistance(queryIndex, index2);
        bool flag_neighbor = true;

        // check for interference
        for (it3 = SR_region.begin(); it3 != SR_region.end(); it3++) {
            unsigned int const index3 = (*it3);

            float const distance13 = sparseMatrix._computeDistance(queryIndex, index3);
            if (distance13 < distance12) {
                float const distance23 = sparseMatrix._computeDistance(index2, index3);
                if (distance23 < distance12) {
                    flag_neighbor = false;
                    break;
                }
            }
        }

        // add neighbor
        if (flag_neighbor) {
            RNG_neighbors[queryIndex].insert(index2);
            RNG_neighbors[index2].insert(queryIndex);
        }
    }
}

// remove RNG links invalidated by the addition of Q
void removeInvalidatedLinks(unsigned int const queryIndex, tsl::sparse_set<unsigned int> const& SR_region, SparseMatrix& sparseMatrix,
                            std::vector<tsl::sparse_set<unsigned int>>& RNG_neighbors) {
    LinksMap checkedLinks;  // don't want to repeat interference checks

    // consider all within region
    tsl::sparse_set<unsigned int>::const_iterator it1, it2;
    for (it1 = SR_region.begin(); it1 != SR_region.end(); it1++) {
        unsigned int const index1 = (*it1);
        float const distance13 = sparseMatrix._computeDistance(index1, queryIndex);

        // consider all neighbors of index1
        tsl::sparse_set<unsigned int> const index1_neighbors = RNG_neighbors[index1];
        for (it2 = index1_neighbors.begin(); it2 != index1_neighbors.end(); it2++) {
            unsigned int const index2 = (*it2);
            if (checkedLinks.checkAndAddLink(index1, index2)) continue;

            // check for query interference
            float const distance12 = sparseMatrix._computeDistance(index1, index2);
            if (distance13 < distance12) {
                float const distance23 = sparseMatrix._computeDistance(index2, queryIndex);
                if (distance23 < distance12) {
                    RNG_neighbors[index1].erase(index2);
                    RNG_neighbors[index2].erase(index1);
                }
            }
        }
    }

    return;
}







void printSet(tsl::sparse_set<unsigned int> const& set) {
    printf("{");
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}

// incremental build of approximate rng by hacid
void approximateRNG(float* const dataPointer, unsigned int const datasetSize, unsigned int const dimension, float epsilon, int numThreads,
                    std::vector<tsl::sparse_set<unsigned int>>& RNG_neighbors, unsigned long long int& distanceComputations, float& time) {
    std::chrono::high_resolution_clock::time_point tStart = std::chrono::high_resolution_clock::now();
    SparseMatrix sparseMatrix(dataPointer, datasetSize, dimension);
    numThreads = std::min(numThreads, omp_get_max_threads());

    RNG_neighbors.clear();
    RNG_neighbors.resize(datasetSize);

    int statsBatchSize = 100;
    printf("%-10s%-10s%-10s\n", "Q", " ave.Dist", " ave.Time");
    float averageDistance = 0;
    std::chrono::duration<double> averageTime;
    unsigned long long int dStart_query = sparseMatrix.get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_query = std::chrono::high_resolution_clock::now();

    unsigned int queryIndex = 0;
    for (queryIndex = 1; queryIndex < datasetSize; queryIndex++) {  // start on 1 since 0 doesn't give edges
        unsigned int const datasetSize = queryIndex;                // easier than a counter

        if (datasetSize == 0) continue;

        // find the nearest neighbor in the graph thus far
        unsigned int nearestNeighbor = 0;
        float nearestDistance = 0;
        nearest_neighbor(queryIndex, datasetSize, sparseMatrix, numThreads, nearestNeighbor, nearestDistance);

        // find furthest link distance of nearestNeighbor
        float const furthestDistance = furthest_neighbor(nearestNeighbor, sparseMatrix, RNG_neighbors);

        // calculate spherical region size
        float SR = (nearestDistance + furthestDistance) * (1 + epsilon);

        // collect all indices within SR
        tsl::sparse_set<unsigned int> SR_region{};
        collectSRIndices(queryIndex, datasetSize, SR, sparseMatrix, numThreads, SR_region);

        // naive RNG neighbors of queryIndex on SR_region
        getRNGNeighbors(queryIndex, SR_region, sparseMatrix, RNG_neighbors);

        // naive invalidation of links by indices in SR_region
        removeInvalidatedLinks(queryIndex, SR_region, sparseMatrix, RNG_neighbors);

        // print stats
        if ((queryIndex + 1) % statsBatchSize == 0) {
            std::chrono::high_resolution_clock::time_point tEnd_query = std::chrono::high_resolution_clock::now();
            unsigned long long int dEnd_query = sparseMatrix.get_distanceComputationCount();
            averageTime = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd_query - tStart_query) / statsBatchSize;
            averageDistance = ((float)(dEnd_query - dStart_query)) / statsBatchSize;
            printf("%-10u%-10.2f%-10.3f\n", queryIndex + 1, averageDistance, averageTime.count() * 1000);
            dStart_query = sparseMatrix.get_distanceComputationCount();
            tStart_query = std::chrono::high_resolution_clock::now();
        }
    }

    std::chrono::high_resolution_clock::time_point tEnd = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();
    distanceComputations = sparseMatrix.get_distanceComputationCount();
    return;
}

// naive nearest neighbor search
void nearest_neighbor1(unsigned int const queryIndex, unsigned int const datasetSize, SearchSparseMatrix& sparseMatrix, int numThreads,
                       unsigned int& nearestNeighbor, float& nearestDistance) {
    nearestDistance = HUGE_VAL;

#pragma omp parallel num_threads(numThreads)
    {
        unsigned int thread_nn = 0;
        float thread_nd = HUGE_VAL;
        SearchSparseMatrixLocal localSparseMatrix = sparseMatrix.getLocalCopy();

#pragma omp for schedule(static)
        for (unsigned int candidateIndex = 0; candidateIndex < datasetSize; candidateIndex++) {
            float const distance = localSparseMatrix._computeDistance(queryIndex, candidateIndex);
            if (distance < thread_nd) {
                thread_nd = distance;
                thread_nn = candidateIndex;
            }
        }

#pragma omp critical(updateNN)
        {
            if (thread_nd < nearestDistance) {
                nearestDistance = thread_nd;
                nearestNeighbor = thread_nn;
            }
            sparseMatrix.updateSparseMatrixWithLocal(localSparseMatrix);
        }
    }

    return;
}

// finding the furthest distance to a neighbor
float furthest_neighbor1(unsigned int const index1, SearchSparseMatrix& sparseMatrix,
                         std::vector<tsl::sparse_set<unsigned int>> const& RNG_neighbors) {
    float furthestDistance = 0;

    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = RNG_neighbors[index1].begin(); it1 != RNG_neighbors[index1].end(); it1++) {
        unsigned int const index2 = (*it1);
        float const distance = sparseMatrix._computeDistance(index1, index2);
        if (distance > furthestDistance) {
            furthestDistance = distance;
        }
    }

    return furthestDistance;
}

// naive nearest neighbor search
void collectSRIndices1(unsigned int const queryIndex, unsigned int const datasetSize, float const SR, SearchSparseMatrix& sparseMatrix,
                       int numThreads, tsl::sparse_set<unsigned int>& SR_region) {
    SR_region.clear();

#pragma omp parallel num_threads(numThreads)
    {
        tsl::sparse_set<unsigned int> thread_SR{};
        SearchSparseMatrixLocal localSparseMatrix = sparseMatrix.getLocalCopy();

#pragma omp for schedule(static)
        for (unsigned int candidateIndex = 0; candidateIndex < datasetSize; candidateIndex++) {
            float const distance = localSparseMatrix._computeDistance(queryIndex, candidateIndex);
            if (distance <= SR) {
                thread_SR.insert(candidateIndex);
            }
        }

#pragma omp critical(updateSR)
        { SR_region.insert(thread_SR.begin(), thread_SR.end()); 
            sparseMatrix.updateSparseMatrixWithLocal(localSparseMatrix);
            }
    }

    return;
}

// get RNG neighbors of query using indices within SR_region
void getRNGNeighbors1(unsigned int const queryIndex, tsl::sparse_set<unsigned int> const& SR_region, SearchSparseMatrix& sparseMatrix,
                      std::vector<tsl::sparse_set<unsigned int>> const& RNG_neighbors, tsl::sparse_set<unsigned int>& searchNeighbors) {
    searchNeighbors.clear();

    // consider all within region as interfering
    tsl::sparse_set<unsigned int>::const_iterator it2, it3;
    for (it2 = SR_region.begin(); it2 != SR_region.end(); it2++) {
        unsigned int const index2 = (*it2);
        float const distance12 = sparseMatrix._computeDistance(queryIndex, index2);
        bool flag_neighbor = true;

        // check for interference
        for (it3 = SR_region.begin(); it3 != SR_region.end(); it3++) {
            unsigned int const index3 = (*it3);

            float const distance13 = sparseMatrix._computeDistance(queryIndex, index3);
            if (distance13 < distance12) {
                float const distance23 = sparseMatrix._computeDistance(index2, index3);
                if (distance23 < distance12) {
                    flag_neighbor = false;
                    break;
                }
            }
        }

        // add neighbor
        if (flag_neighbor) {
            searchNeighbors.insert(index2);
        }
    }
    return;
}

// incremental build of approximate rng by hacid
void search(unsigned int const dimension, float* const dataPointer, unsigned int const datasetSize, float* const testPointer,
            unsigned int const testsetSize, float epsilon, int numThreads, std::vector<tsl::sparse_set<unsigned int>>& RNG_neighbors,
            std::vector<tsl::sparse_set<unsigned int>>& searchNeighbors, float& distanceComputations, float& time) {
    std::chrono::high_resolution_clock::time_point tStart = std::chrono::high_resolution_clock::now();
    std::shared_ptr<SparseMatrix> sparseMatrix1 = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension);
    SearchSparseMatrix sparseMatrix = SearchSparseMatrix(testPointer, testsetSize, dimension, sparseMatrix1);
    numThreads = std::min(numThreads, omp_get_max_threads());

    searchNeighbors.clear();
    searchNeighbors.resize(testsetSize);

    int statsBatchSize = 100;
    printf("%-10s%-10s%-10s\n", "Q", " ave.Dist", " ave.Time");
    float averageDistance = 0;
    std::chrono::duration<double> averageTime;
    unsigned long long int dStart_query = sparseMatrix.get_distanceComputationCount();
    std::chrono::high_resolution_clock::time_point tStart_query = std::chrono::high_resolution_clock::now();

    unsigned int queryIndex = 0;
    for (queryIndex = 0; queryIndex < testsetSize; queryIndex++) {  // start on 1 since 0 doesn't give edges
        unsigned int const index1 = datasetSize + queryIndex;

        // find the nearest neighbor in the graph thus far
        unsigned int nearestNeighbor = 0;
        float nearestDistance = 0;
        nearest_neighbor1(index1, datasetSize, sparseMatrix, numThreads, nearestNeighbor, nearestDistance);

        // find furthest link distance of nearestNeighbor
        float const furthestDistance = furthest_neighbor1(nearestNeighbor, sparseMatrix, RNG_neighbors);

        // calculate spherical region size
        float SR = (nearestDistance + furthestDistance) * (1 + epsilon);

        // collect all indices within SR
        tsl::sparse_set<unsigned int> SR_region{};
        collectSRIndices1(index1, datasetSize, SR, sparseMatrix, numThreads, SR_region);

        // naive RNG neighbors of queryIndex on SR_region
        searchNeighbors[queryIndex] = tsl::sparse_set<unsigned int>{};
        getRNGNeighbors1(index1, SR_region, sparseMatrix, RNG_neighbors, searchNeighbors[queryIndex]);
    }

    std::chrono::high_resolution_clock::time_point tEnd = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() * 1000 / testsetSize;
    distanceComputations = sparseMatrix.get_distanceComputationCount() / testsetSize;
    return;
}

}  // namespace hacid
