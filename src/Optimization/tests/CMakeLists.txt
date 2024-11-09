add_executable(test_optimization test_optimization.cpp)

# linking optimization
target_include_directories(test_optimization  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/Optimization/Optimization/)
target_link_libraries(test_optimization PRIVATE optimization)

# link to dataset creator, header only
target_include_directories(test_optimization  PRIVATE 
        ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG/)

# linking boost unit testing
target_include_directories(test_optimization        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_optimization PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

# add to ctest
add_test(NAME test_optimization COMMAND test_optimization)

