set(HEADERS pivot-layer-hierarchy.hpp pivot-layer.hpp index-pair-distance-struct.hpp)
set(SOURCES pivot-layer-hierarchy.cpp 
            pivot-layer.cpp
            umax.cpp
            vmax.cpp)

add_library(pivot_layer_hierarchy ${HEADERS} ${SOURCES})

# link to sparse matrix 
target_include_directories(pivot_layer_hierarchy PUBLIC ${PROJECT_SOURCE_DIR}/src/SparseMatrix/SparseMatrix)
target_link_libraries(pivot_layer_hierarchy PUBLIC sparse_matrix)

# link to tsl 
target_include_directories(pivot_layer_hierarchy PUBLIC ${PROJECT_SOURCE_DIR}/extern/sparse-map/include)
