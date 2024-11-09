#define BOOST_TEST_DYN_LINK  // deals with no main() issue

#include "Test_IncrementalBuild.hpp"

BOOST_AUTO_TEST_SUITE(incremental_build)

BOOST_AUTO_TEST_CASE(coarseNeighbors) { Test_IncrementalBuild::test_coarseNeighbors(); }

BOOST_AUTO_TEST_CASE(vmax) { Test_IncrementalBuild::test_vmax(); }

BOOST_AUTO_TEST_CASE(pivot_selection) { Test_IncrementalBuild::test_pivotSelection(); }

BOOST_AUTO_TEST_CASE(family_tree) { Test_IncrementalBuild::test_familyTree(); }

BOOST_AUTO_TEST_CASE(dmax) { Test_IncrementalBuild::test_dmax(); }

BOOST_AUTO_TEST_CASE(umax) { Test_IncrementalBuild::test_umax(); }

BOOST_AUTO_TEST_CASE(build) { Test_IncrementalBuild::test_neighbors(); }

BOOST_AUTO_TEST_SUITE_END()
