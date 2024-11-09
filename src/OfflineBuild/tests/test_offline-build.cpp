#define BOOST_TEST_DYN_LINK  // deals with no main() issue

#include "Test_OfflineBuild.hpp"

BOOST_AUTO_TEST_SUITE(offline_build)

BOOST_AUTO_TEST_CASE(coarse_neighbors) {
    Test_OfflineBuild::test_coarseNeighbors();
}

BOOST_AUTO_TEST_CASE(family_tree) {
    Test_OfflineBuild::test_familyTree();
}

BOOST_AUTO_TEST_CASE(dmax) {
    Test_OfflineBuild::test_dmax();
}

BOOST_AUTO_TEST_CASE(umax) {
    Test_OfflineBuild::test_umax();
}

BOOST_AUTO_TEST_CASE(build) {
    Test_OfflineBuild::test_neighbors();
}



BOOST_AUTO_TEST_SUITE_END()
