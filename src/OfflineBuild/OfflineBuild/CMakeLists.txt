set(HEADERS offline-build.hpp query-struct-offline-build.hpp)
set(SOURCES offline-build.cpp
            stage0.cpp 
            stage1.cpp
            stage2.cpp 
            stage3.cpp 
            stage4.cpp
            stage5.cpp 
            stage6.cpp 
            stage7.cpp)

add_library(offline_build ${HEADERS} ${SOURCES})

# add open mp 
target_link_libraries(offline_build PUBLIC OpenMP::OpenMP_CXX)

# pivot selection
add_subdirectory(PivotSelection/PivotSelection)
target_include_directories(offline_build PUBLIC PivotSelection/PivotSelection)
target_link_libraries(offline_build PUBLIC pivot_selection)
if(IHRNG_BUILD_TESTING)
    add_subdirectory(PivotSelection/tests)
endif()

# link to pivot layer hierarchy
target_include_directories(offline_build PUBLIC ${PROJECT_SOURCE_DIR}/src/PivotLayerHierarchy/PivotLayerHierarchy)
target_link_libraries(offline_build PUBLIC pivot_layer_hierarchy)

# link to Utils
target_include_directories(offline_build PUBLIC ${PROJECT_SOURCE_DIR}/src/Utils/Utils)
target_link_libraries(offline_build PRIVATE memory_usage)

