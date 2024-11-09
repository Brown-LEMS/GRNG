set(HEADERS1 Test-PivotLayer.hpp)
set(SOURCES1 Test-PivotLayer.cpp)
add_executable(test_pivot_layer             test_pivot-layer.cpp            ${HEADERS1} ${SOURCES1})

set(HEADERS2 Test-PivotLayerHierarchy.hpp)
set(SOURCES2 Test-PivotLayerHierarchy.cpp)
add_executable(test_pivot_layer_hierarchy   test_pivot-layer-hierarchy.cpp  ${HEADERS2} ${SOURCES2})

# linking pivot-layer-hierarchy library
target_include_directories(test_pivot_layer             PRIVATE ${PROJECT_SOURCE_DIR}/src/PivotLayerHierarchy/PivotLayerHierarchy)
target_include_directories(test_pivot_layer_hierarchy   PRIVATE ${PROJECT_SOURCE_DIR}/src/PivotLayerHierarchy/PivotLayerHierarchy)
target_link_libraries(test_pivot_layer                  PRIVATE pivot_layer_hierarchy)
target_link_libraries(test_pivot_layer_hierarchy        PRIVATE pivot_layer_hierarchy)

# linking boost unit testing
target_include_directories(test_pivot_layer        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_pivot_layer                  PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})
target_include_directories(test_pivot_layer_hierarchy        PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_pivot_layer_hierarchy        PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

# add to ctest
add_test(NAME test_pivot_layer              COMMAND test_pivot_layer)
add_test(NAME test_pivot_layer_hierarchy    COMMAND test_pivot_layer_hierarchy)

