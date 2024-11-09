#define BOOST_TEST_DYN_LINK  // deals with no main() issue
#include "Test-SparseMatrix.hpp"


BOOST_AUTO_TEST_SUITE(sparse_matrix)

BOOST_AUTO_TEST_CASE(compute_distance) {
    Test_SparseMatrix::test_computeDistance();
}

BOOST_AUTO_TEST_CASE(insert_distance) {
    Test_SparseMatrix::test_insertDistance();
}

BOOST_AUTO_TEST_CASE(get_distance) {
    Test_SparseMatrix::test_getDistance();
}

BOOST_AUTO_TEST_CASE(local_matrix) {
    Test_SparseMatrix::test_localMatrix();
}

BOOST_AUTO_TEST_SUITE_END()