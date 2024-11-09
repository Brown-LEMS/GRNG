set(HEADERS online-search.hpp
            query-struct-online-search.hpp)

set(SOURCES online-search.cpp 
            brute-force.cpp
            stage0.cpp
            stage1.cpp
            stage2.cpp
            stage3.cpp
            stage4.cpp
            stage5.cpp
            stage6.cpp)

add_library(online_search ${HEADERS} ${SOURCES})

# add open mp 
target_link_libraries(online_search PUBLIC OpenMP::OpenMP_CXX)

# link to pivot layer hierarchy
target_include_directories(online_search PUBLIC ${PROJECT_SOURCE_DIR}/src/PivotLayerHierarchy/PivotLayerHierarchy)
target_link_libraries(online_search PUBLIC pivot_layer_hierarchy)

# add search sparse matrix
target_include_directories(online_search PUBLIC ${PROJECT_SOURCE_DIR}/src/SparseMatrix/SparseMatrix)
target_link_libraries(online_search PUBLIC search_sparse_matrix)

# link to Utils
target_include_directories(online_search PUBLIC ${PROJECT_SOURCE_DIR}/src/Utils/Utils)
target_link_libraries(online_search PRIVATE memory_usage)
