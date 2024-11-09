set(HEADERS Test_OfflineBuild.hpp)
set(SOURCES Test_OfflineBuild.cpp)

add_executable(test_offline_build test_offline-build.cpp ${HEADERS} ${SOURCES})

# linking offline-build
target_include_directories(test_offline_build  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/OfflineBuild/OfflineBuild/)
target_link_libraries(test_offline_build  PRIVATE offline_build)

# link to dataset creator, header only
target_include_directories(test_offline_build  PRIVATE 
        ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG/)

# linking boost unit testing
target_include_directories(test_offline_build        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_offline_build PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

# add to ctest
add_test(NAME test_offline_build COMMAND test_offline_build)


