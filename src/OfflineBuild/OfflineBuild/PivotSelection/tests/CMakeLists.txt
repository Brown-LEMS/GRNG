add_executable(test_pivot_selection test_pivot-selection.cpp )

# linking pivot-selection 
target_include_directories(test_pivot_selection  PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/OfflineBuild/OfflineBuild/PivotSelection/PivotSelection/)
target_link_libraries(test_pivot_selection  PRIVATE pivot_selection)

# link to dataset creator, header only
target_include_directories(test_pivot_selection  PRIVATE 
        ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG/)

# linking boost unit testing
target_link_libraries(test_pivot_selection PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

# add to ctest
add_test(NAME test_pivot_selection COMMAND test_pivot_selection)

