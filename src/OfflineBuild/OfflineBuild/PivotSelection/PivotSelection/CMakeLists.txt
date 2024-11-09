set(HEADERS constant-radius-pivot-selection.hpp 
            constant-radius-pivot-selection-struct.hpp
            pivot-selection-struct.hpp)
set(SOURCES constant-radius-pivot-selection.cpp)
add_library(pivot_selection ${HEADERS} ${SOURCES})

# link to sparse matrix 
target_include_directories(pivot_selection PUBLIC ${PROJECT_SOURCE_DIR}/src/SparseMatrix/SparseMatrix)
target_link_libraries(pivot_selection PUBLIC sparse_matrix)

# link to tsl 
target_include_directories(pivot_selection PUBLIC ${PROJECT_SOURCE_DIR}/extern/sparse-map/include)
