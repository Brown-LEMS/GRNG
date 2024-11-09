#define BOOST_TEST_DYN_LINK  // deals with no main() issue

#include "Test_OnlineSearch.hpp"

BOOST_AUTO_TEST_SUITE(online_search)

BOOST_AUTO_TEST_CASE(neighbors) { Test_OnlineSearch::test_neighbors(); }

BOOST_AUTO_TEST_SUITE_END()
