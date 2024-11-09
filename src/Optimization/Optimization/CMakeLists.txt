set(HEADERS optimization.hpp
            hrng-data.hpp
            hrng-results.hpp)
set(SOURCES optimization.cpp)

add_library(optimization ${HEADERS} ${SOURCES})

#link nlopt library for optimization
target_include_directories(optimization PRIVATE ${NLOPT_DIRECTORIES})
target_link_libraries(optimization PRIVATE ${NLOPT})

# link to pivot layer hierarchy
target_include_directories(optimization PUBLIC ${PROJECT_SOURCE_DIR}/src/PivotLayerHierarchy/PivotLayerHierarchy)
target_link_libraries(optimization PUBLIC pivot_layer_hierarchy)

# link to IncrementalBuild
target_include_directories(optimization PUBLIC ${PROJECT_SOURCE_DIR}/src/IncrementalBuild/IncrementalBuild)
target_link_libraries(optimization PUBLIC incremental_build)

# link to OnlineSearch
target_include_directories(optimization PUBLIC ${PROJECT_SOURCE_DIR}/src/OnlineSearch/OnlineSearch)
target_link_libraries(optimization PUBLIC online_search)

# link to ResultsSaver
target_include_directories(optimization PUBLIC ${PROJECT_SOURCE_DIR}/src/Utils/Utils)
target_link_libraries(optimization PUBLIC results_saver)
