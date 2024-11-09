#ifndef Test_PivotLayerHierarchy_hpp
#define Test_PivotLayerHierarchy_hpp

#define BOOST_TEST_MODULE PivotLayerHierarchyTest
#include <boost/test/unit_test.hpp>

#include "pivot-layer-hierarchy.hpp"

namespace Test_PivotLayerHierarchy {
void test_constructors();
void test_sparseMatrix();
void test_ancestryDotCom();
void test_umax();
void test_umaxDescendants();
void test_umaxMaximumOfDescendants();
void test_copyConstructor();
};  // namespace Test_PivotLayerHierarchy

#endif  // Test_PivotLayerHierarchy_hpp