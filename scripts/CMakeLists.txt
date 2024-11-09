# Simple build 
add_executable(build build.cpp)
target_link_libraries(build PRIVATE hRNG)
target_include_directories(build PRIVATE ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG)
install(TARGETS build)

# Build and Search
add_executable(build_search build-search.cpp)
target_link_libraries(build_search PRIVATE hRNG)
target_include_directories(build_search PRIVATE ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG)
install(TARGETS build_search)

# Build and Search
add_executable(testResults testResults.cpp)
target_link_libraries(testResults PRIVATE hRNG)
target_include_directories(testResults PRIVATE ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG)
install(TARGETS testResults)

# Pivot Selection
add_executable(selectPivots selectPivots.cpp)
target_link_libraries(selectPivots PRIVATE hRNG)
target_include_directories(selectPivots PRIVATE ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG)
install(TARGETS selectPivots)

# Optimization
if(HRNG_OPT)
    add_executable(optimize optimize.cpp)
    target_include_directories(optimize PRIVATE ${PROJECT_SOURCE_DIR}/include/HierarchicalRNG)
    target_include_directories(optimize PRIVATE ${PROJECT_SOURCE_DIR}/src/Optimization/Optimization/)
    target_link_libraries(optimize PRIVATE optimization)
    install(TARGETS optimize)
endif()
