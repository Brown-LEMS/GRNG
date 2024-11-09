#include "Test-SparseMatrix.hpp"

#include <iostream>

/*  test_computeDistance()
      - create a float pointer with 3,4,5 relationship
      - initialize sparse matrix
      - compute distances
      - check if distance is correct
*/
void Test_SparseMatrix::test_computeDistance() {
    bool success = true;

    float data[4] = {6, 8, 3, 4};
    float* dataPointer = data;
    SparseMatrix SM(dataPointer, 2, 2);

    // ensure correct computation
    if (SM._computeDistance(0, 1) != 5.0f) {
        BOOST_FAIL("Incorrect Euclidean Distance");
    } else if (SM.get_distanceComputationCount() != 1) {
        BOOST_FAIL("Incorrect Distance Count");
    }

    return;
}

/*  test_insertDistance()
      - create a float pointer
      - initialize sparse matrix
      - store arbitrary distance between two points
      - retrieve using getDistanceIfAvailable
      - check both sides, (0,1) and (1,0)
      - ensure both return true for distance present, and distance is the same
*/
void Test_SparseMatrix::test_insertDistance() {
    bool success = true;

    float data[10] = {6, 8, 3, 4, 3, 6, 4, 3, 4, 3};
    float* dataPointer = data;
    SparseMatrix SM(dataPointer, 5, 2);

    float initial = 5.0f;
    SM._insertDistanceIntoMatrix(0, 1, initial);
    std::pair<bool, float> retreived1 = SM.getDistanceIfAvailable(0, 1);
    std::pair<bool, float> retreived2 = SM.getDistanceIfAvailable(1, 0);

    // check distances are there
    if (!retreived1.first || !retreived2.first) {
        BOOST_FAIL("retrieve distance failed");
    }

    // check distances are equal
    if ((retreived1.second != initial) || (retreived2.second != initial)) {
        BOOST_FAIL("incorrect retrieved distance value");
    }

    return;
}

/*  test_getDistance()
      - create a float pointer
      - initialize sparse matrix
      - store arbitrary distance between two points
      - retrieve using getDistanceIfAvailable
      - check both sides, (0,1) and (1,0)
      - ensure both return true for distance present, and distance is the same
*/
void Test_SparseMatrix::test_getDistance() {
    bool success = true;

    float data[10] = {6, 8, 3, 4, 3, 6, 4, 3, 4, 3};
    float* dataPointer = data;
    SparseMatrix SM(dataPointer, 5, 2);

    // ensure distance computatoin count starts at zero
    if (SM.get_distanceComputationCount() != 0) {
        BOOST_FAIL("sparsematrix initialized with nonzero distance count");
    }

    // check correct computation
    float distance = SM.getDistance(0, 1);
    if (distance != (float)5.0f) {
        BOOST_FAIL("euclidean distance incorrect");
    } else if (SM.get_distanceComputationCount() != 1) {
        BOOST_FAIL("distance count incorrect");
    }

    // get distance again
    float distance1 = SM.getDistance(0, 1);
    if (distance != distance1) {
        BOOST_FAIL("incorrect distance fetch");
    } else if (SM.get_distanceComputationCount() != 1) {
        BOOST_FAIL("distance count incorrect");
    }

    // get distance other side
    float distance2 = SM.getDistance(1, 0);
    if (distance != distance2) {
        BOOST_FAIL("distances not symmetric");
    } else if (SM.get_distanceComputationCount() != 1) {
        BOOST_FAIL("distance count incorrect");
    }

    // get new distance, but same index
    // shouldn't compute distance, just return 0
    float const distance3 = SM.getDistance(2, 2);
    if (distance3 != 0.0f) {
        BOOST_FAIL("non-zero same index distance");
    } else if (SM.get_distanceComputationCount() != 1) {
        BOOST_FAIL("same-index distance incorrectly computed");
    }

    return;
}

void Test_SparseMatrix::test_localMatrix() {
    float data[10] = {6, 8, 3, 4, 3, 6, 4, 3, 4, 3};
    float* dataPointer = data;
    SparseMatrix SM(dataPointer, 5, 2);

    // get one distance
    float const distance1a = SM.getDistance(0, 1);

    // get local copy of sparse matrix
    SparseMatrixLocal localSM = SM.getLocalCopy();
    std::pair<bool, float> const distance1b = SM.getDistanceIfAvailable_LocalAndGlobal(0, 1, localSM);
    if (!distance1b.first) {
        BOOST_FAIL("distance not available by local");
    }
    if (distance1b.second != distance1a) {
        BOOST_FAIL("Distance not fetched from global");
    }

    // compute distance from local
    float const distance3b = SM.getDistance_LocalAndGlobal(3, 1, localSM);
    std::pair<bool, float> const distance3a = SM.getDistanceIfAvailable(3, 1);  // shouldnt be available
    if (distance3a.first) {
        BOOST_FAIL("local distance computation stored in global");
    }

    // merge
    SM.updateSparseMatrixWithLocal(localSM);
    std::pair<bool, float> const distance3c = SM.getDistanceIfAvailable(3, 1);  // shouldnt be available
    if (distance3c.first == false) {
        BOOST_FAIL("global matrix not updated with local!");
    }

    return;
}