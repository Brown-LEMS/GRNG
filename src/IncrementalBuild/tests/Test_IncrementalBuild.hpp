#ifndef Test_IncrementalBuild_hpp
#define Test_IncrementalBuild_hpp

#define BOOST_TEST_MODULE IncrementalBuildTest


#include <boost/test/unit_test.hpp>

#include "incremental-build.hpp"

namespace Test_IncrementalBuild {

void test_pivotSelection();
void test_familyTree();
void test_dmax();
void test_umax();
void test_neighbors();
void test_coarseNeighbors();
void test_vmax();
}  // namespace Test_IncrementalBuild

#endif  // Test_IncrementalBuild_hpp