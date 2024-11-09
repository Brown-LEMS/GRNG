set(HEADERS Test_IncrementalBuild.hpp)
set(SOURCES Test_IncrementalBuild.cpp)

add_executable(test_incremental_build test_incremental-build.cpp ${HEADERS} ${SOURCES})

# linking incremental-build
target_include_directories(test_incremental_build  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/IncrementalBuild/IncrementalBuild/)
target_link_libraries(test_incremental_build  PRIVATE incremental_build)

# linking offline-build for verification
target_include_directories(test_incremental_build  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/OfflineBuild/OfflineBuild/)
target_link_libraries(test_incremental_build  PRIVATE offline_build)

# link to dataset creator, header only
target_include_directories(test_incremental_build  PRIVATE 
        ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG/)

# linking boost unit testing
target_include_directories(test_incremental_build        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_incremental_build PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

# add to ctest
add_test(NAME test_incremental_build COMMAND test_incremental_build)

