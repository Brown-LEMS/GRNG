#ifndef Test_SparseMatrix_hpp
#define Test_SparseMatrix_hpp

#define BOOST_TEST_MODULE SparseMatrixTest
#include <boost/test/unit_test.hpp>

#include "sparse-matrix.hpp"


namespace Test_SparseMatrix {

// compute distance test
void test_computeDistance();

// insert distance test
void test_insertDistance();

// get distance test
void test_getDistance();


void test_localMatrix();

};  // namespace Test_SparseMatrix

#endif  // Test_SparseMatrix_hpp