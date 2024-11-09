# library for SparseMatrix
set(HEADERS sparse-matrix.hpp
            sparse-matrix-local.hpp)
set(SOURCES sparse-matrix.cpp
            sparse-matrix-local.cpp)
add_library(sparse_matrix ${HEADERS} ${SOURCES})

# hash maps
target_include_directories(sparse_matrix PUBLIC ${PROJECT_SOURCE_DIR}/extern/flat-hash-map)
target_include_directories(sparse_matrix PUBLIC ${PROJECT_SOURCE_DIR}/extern/sparse-map/include)

# library for SearchSparseMatrix
set(HEADERS search-sparse-matrix.hpp
            search-sparse-matrix-local.hpp)
set(SOURCES search-sparse-matrix.cpp
            search-sparse-matrix-local.cpp)
add_library(search_sparse_matrix ${HEADERS} ${SOURCES})

# link sparse matrix
target_link_libraries(search_sparse_matrix PUBLIC sparse_matrix)