set(HEADERS Test_OnlineSearch.hpp)
set(SOURCES Test_OnlineSearch.cpp)

add_executable(test_online_search test_online-search.cpp ${HEADERS} ${SOURCES})

# linking online-search
target_include_directories(test_online_search  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/OnlineSearch/OnlineSearch/)
target_link_libraries(test_online_search PRIVATE online_search)

# linking offline-build for verification
target_include_directories(test_online_search  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/OfflineBuild/OfflineBuild/)
target_link_libraries(test_online_search  PRIVATE offline_build)

# linking incremental-build for verification
target_include_directories(test_online_search  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/IncrementalBuild/IncrementalBuild/)
target_link_libraries(test_online_search  PRIVATE incremental_build)

# link to dataset creator, header only
target_include_directories(test_online_search  PRIVATE 
        ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG/)

# linking boost unit testing
target_include_directories(test_online_search        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_online_search PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

# add to ctest
add_test(NAME test_online_search COMMAND test_online_search)

