#define BOOST_TEST_DYN_LINK  // deals with no main() issue

#include "Test-PivotLayerHierarchy.hpp"

BOOST_AUTO_TEST_SUITE(pivot_layer_hierarchy)

BOOST_AUTO_TEST_CASE(constructors) {
    Test_PivotLayerHierarchy::test_constructors();
}

BOOST_AUTO_TEST_CASE(copyConstructor) {
    Test_PivotLayerHierarchy::test_copyConstructor();
}

BOOST_AUTO_TEST_CASE(sparseMatrix) {
    Test_PivotLayerHierarchy::test_sparseMatrix();
}

BOOST_AUTO_TEST_CASE(ancestryDotCom) {
    Test_PivotLayerHierarchy::test_ancestryDotCom();
}

BOOST_AUTO_TEST_CASE(umax) {
    Test_PivotLayerHierarchy::test_umax();
}

BOOST_AUTO_TEST_SUITE_END()
