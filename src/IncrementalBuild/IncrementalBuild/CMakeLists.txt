set(HEADERS incremental-build.hpp
            query-struct-incremental-build.hpp)
set(SOURCES incremental-build.cpp 
            brute-force.cpp
            stage0.cpp
            stage1.cpp
            stage2.cpp
            stage3.cpp
            stage4.cpp
            stage5.cpp
            stage6.cpp
            stage7.cpp
            stage8.cpp
            stage9.cpp)

add_library(incremental_build ${HEADERS} ${SOURCES})

# add open mp 
target_link_libraries(incremental_build PUBLIC OpenMP::OpenMP_CXX)

# link to pivot layer hierarchy
target_include_directories(incremental_build PUBLIC ${PROJECT_SOURCE_DIR}/src/PivotLayerHierarchy/PivotLayerHierarchy)
target_link_libraries(incremental_build PUBLIC pivot_layer_hierarchy)

# link to Utils
target_include_directories(incremental_build PUBLIC ${PROJECT_SOURCE_DIR}/src/Utils/Utils)
target_link_libraries(offline_build PRIVATE memory_usage)
