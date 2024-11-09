#ifndef Test_PivotLayer_hpp
#define Test_PivotLayer_hpp

#define BOOST_TEST_MODULE PivotLayerTest
#include <boost/test/unit_test.hpp>

#include "pivot-layer.hpp"

namespace Test_PivotLayer {
void test_constructors();
void test_addPivotToLayer();
void test_initializePivotsInLayer();
void test_sharedPointers();
void test_sparseMatrix();
void test_copyConstructor();
void test_lune();
void test_bruteForceConstruction();
}  // namespace Test_PivotLayer

#endif  // Test_PivotLayer_hpp