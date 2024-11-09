#ifndef Test_OfflineBuild_hpp
#define Test_OfflineBuild_hpp

#define BOOST_TEST_MODULE OfflineBuildTest

#include <boost/test/unit_test.hpp>

#include "offline-build.hpp"

namespace Test_OfflineBuild {

void test_familyTree();
void test_dmax();
void test_umax();
void test_neighbors();
void test_coarseNeighbors();


}  // namespace Test_OfflineBuild

#endif  // Test_OfflineBuild_hpp