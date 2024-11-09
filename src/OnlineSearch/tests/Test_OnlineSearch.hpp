#ifndef Test_OnlineSearch_hpp
#define Test_OnlineSearch_hpp

#define BOOST_TEST_MODULE OnlineSearchTest

#include <boost/test/unit_test.hpp>

#include "online-search.hpp"

namespace Test_OnlineSearch {



void getRNGNeighbors(float* dataPointer, unsigned int const datasetSize, unsigned int const dimension, float* testPointer,
                     unsigned int const testsetSize, tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& RNG);
void test_neighbors();



inline std::vector<float> getEffectiveRadiusVector(std::vector<float> vec) {
    std::vector<float> new_vec = vec;
    for (int i = 0; i < (int)vec.size(); i++) {
        for (int j = i + 1; j < (int)vec.size(); j++) {
            new_vec[i] += (float)vec[j];
        }
    }
    return new_vec;
}

inline void printSet(tsl::sparse_set<unsigned int> const& set) {
    std::cout << "{";
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        std::cout << (*it1) << ",";
    }
    std::cout << "}, \n";
}

}  // namespace Test_OnlineSearch

#endif  // Test_OnlineSearch_hpp