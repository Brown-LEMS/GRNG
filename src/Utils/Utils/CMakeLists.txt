# rank-ordered list
add_library(rank_ordered_list rank-ordered-list.hpp index-distance-struct.hpp)
set_target_properties(rank_ordered_list PROPERTIES LINKER_LANGUAGE CXX)

# links-map
add_library(links_map links-map.hpp)
set_target_properties(links_map PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(links_map PUBLIC ${PROJECT_SOURCE_DIR}/extern/sparse-map/include)

# memory-usage
add_library(memory_usage memory-usage.h memory-usage.cpp)
set_target_properties(memory_usage PROPERTIES LINKER_LANGUAGE CXX)

# results saver
add_library(results_saver results-saver.hpp)
set_target_properties(results_saver PROPERTIES LINKER_LANGUAGE CXX)