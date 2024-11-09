add_library(hRNG OBJECT hierarchical-rng.cpp)

# need to retrieve header file from include/ directory
target_include_directories(hRNG PRIVATE "${PROJECT_SOURCE_DIR}/include/HierarchicalRNG/")

# Utils subdirectory 
add_subdirectory(Utils/Utils)

# include SparseMatrix subdirectory to build library first
add_subdirectory(SparseMatrix/SparseMatrix)
if(HRNG_BUILD_TESTING)
    add_subdirectory(SparseMatrix/tests)
endif()

# include PivotLayerHierarchy subdirectory
add_subdirectory(PivotLayerHierarchy/PivotLayerHierarchy)
target_include_directories(hRNG PUBLIC PivotLayerHierarchy/PivotLayerHierarchy)
target_link_libraries(hRNG PUBLIC pivot_layer_hierarchy)
if(HRNG_BUILD_TESTING)
    add_subdirectory(PivotLayerHierarchy/tests)
endif()

# include OfflineBuild subdirectory
add_subdirectory(OfflineBuild/OfflineBuild)
target_include_directories(hRNG PUBLIC OfflineBuild/OfflineBuild)
target_link_libraries(hRNG PUBLIC offline_build)
if(HRNG_BUILD_TESTING)
    add_subdirectory(OfflineBuild/tests)
endif()

# include IncrementalBuild subdirectory
add_subdirectory(IncrementalBuild/IncrementalBuild)
target_include_directories(hRNG PUBLIC IncrementalBuild/IncrementalBuild)
target_link_libraries(hRNG PUBLIC incremental_build)
if(HRNG_BUILD_TESTING)
    add_subdirectory(IncrementalBuild/tests)
endif()

# include OnlineSearch subdirectory
add_subdirectory(OnlineSearch/OnlineSearch)
target_include_directories(hRNG PUBLIC OnlineSearch/OnlineSearch)
target_link_libraries(hRNG PUBLIC online_search)
if(HRNG_BUILD_TESTING)
    add_subdirectory(OnlineSearch/tests)
endif()

# include Optimization
if(HRNG_OPT)
    add_subdirectory(Optimization/Optimization)
    target_include_directories(hRNG PUBLIC Optimization/Optimization)
    target_link_libraries(hRNG PUBLIC optimization)
    if(HRNG_BUILD_TESTING)
        add_subdirectory(Optimization/tests)
    endif()
endif()

