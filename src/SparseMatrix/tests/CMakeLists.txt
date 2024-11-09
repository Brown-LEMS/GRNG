set(HEADERS Test-SparseMatrix.hpp)
set(SOURCES Test-SparseMatrix.cpp)
add_executable(test_sparse_matrix test-sparse-matrix.cpp ${HEADERS} ${SOURCES})

target_include_directories(test_sparse_matrix PRIVATE ${PROJECT_SOURCE_DIR}/src/SparseMatrix/SparseMatrix)
target_link_libraries(test_sparse_matrix PRIVATE sparse_matrix)

# linking boost unit testing
target_include_directories(test_sparse_matrix        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_sparse_matrix        PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

add_test(NAME test_sparse_matrix COMMAND test_sparse_matrix)