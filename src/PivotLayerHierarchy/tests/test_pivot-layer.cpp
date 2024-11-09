#define BOOST_TEST_DYN_LINK  // deals with no main() issue

#include "Test-PivotLayer.hpp"

BOOST_AUTO_TEST_SUITE(pivot_layer)

BOOST_AUTO_TEST_CASE(constructors) {
    Test_PivotLayer::test_constructors();
}

BOOST_AUTO_TEST_CASE(addPivotToLayer) {
    Test_PivotLayer::test_addPivotToLayer();
}

BOOST_AUTO_TEST_CASE(initializePivotsInLayer) {
    Test_PivotLayer::test_initializePivotsInLayer();
}

BOOST_AUTO_TEST_CASE(sharedPointers) {
    Test_PivotLayer::test_sharedPointers();
}

BOOST_AUTO_TEST_CASE(sparseMatrix) {
    Test_PivotLayer::test_sparseMatrix();
}

BOOST_AUTO_TEST_CASE(copyConstructors) {
    Test_PivotLayer::test_copyConstructor();
}

BOOST_AUTO_TEST_CASE(lune) {
    Test_PivotLayer::test_lune();
}

BOOST_AUTO_TEST_CASE(bruteForceConstruction) {
    Test_PivotLayer::test_bruteForceConstruction();
}

BOOST_AUTO_TEST_SUITE_END()
