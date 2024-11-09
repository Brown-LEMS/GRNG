cmake_minimum_required(VERSION 3.15)
project(hRNG 
        VERSION 1.0
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)  
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# User defined options
option(CMAKE_BUILD_TYPE "type of build [Release/Debug/RelWithDebInfo]" Release)
option(BUILD_TESTING "unit testing" OFF)
option(HRNG_BUILD_TESTING "unit testing for ihRNG specifically (if external usage)" OFF)
option(HRNG_OPT "enable use/build of optimization library" OFF)

# openmp
find_package(OpenMP)

# unit testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message("Build Testing: ${BUILD_TESTING}")
    if (BUILD_TESTING) # ensure main project is ihRNG
        set(HRNG_BUILD_TESTING ON)
    endif()
endif()
if(HRNG_BUILD_TESTING)
    set(HRNG_BUILD_TESTING ON)
    include(CTest)

    # get boost testing library
    find_package(Boost COMPONENTS unit_test_framework REQUIRED)
endif()

# hRNG optimization
if (HRNG_OPT)
    #SET(NLOPT_HOME "/home/colef8/nlopt")
    SET(NLOPT_HOME "/users/cfoste18/data/cfoste18/hRNG/hRNG_Optimization/include/nlopt/install")
    SET(NLOPT_DIRECTORIES "${NLOPT_HOME}/include/")
    find_library(NLOPT nlopt PATHS "${NLOPT_HOME}/lib/")
endif()
MESSAGE("NLOPT: ${NLOPT}")

# Hierarchical RNG
add_subdirectory(src)

# set CMAKE_INSTALL_PREFIX as default location of `make install` command
# creates a bin/ folder within the location
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# add scripts
add_subdirectory(scripts)