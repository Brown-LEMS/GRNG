# Generalized Relative Neighborhood Graph (GRNG)

### Overview:
The Generalized RNG (GRNG) is an incremental, hierarchical approach to constructing the Relative Neighborhood Graph (RNG).

Each layer of the hierarchy contains a graph that is a generalization of the RNG, where each node is a "pivot" in that layer with an associated radius. The bottom, RNG layer includes all exemplar has a radius of 0. 
- The radius of a pivot determines the size of its pivot domain (hypersphere around the pivot), and pivots in lower layers within its domain are considered children of that pivot. 
- This uses an incremental approach to building the index. 

## Requirements and Dependencies
- Tested on gcc 8.3.0 and CMake 3.15
- Uses ![TSL](https://github.com/Tessil/sparse-map) as a hash map. 
    - `cd extern/ && git clone https://github.com/Tessil/sparse-map.git sparse-map/`
- Uses BOOST library for unit testing and serialization
    - module load boost/1.68
- Uses NLOPT for optimization

## Build and Run 
cmake -S . -B build/
cd build ; make install
cd bin ; ./main

## Unit Testing
Boost is used for unit testing. 
- unit testing turned off by default, but recommended to run testing once to ensure everything is working properly
- to turn on unit testing, change `option(BUILD_TESTING "unit testing" OFF)` to `ON` within CMakeLists.txt before compiling
    - or use `ccmake .` to change cmake env variable after compiling
- to run unit testing, compile and in build folder run:
    - `ctest` or `ctest --output-on-failure` or `ctest --extra-verbose`

## Optimization
![NLOPT](https://nlopt.readthedocs.io/en/latest/NLopt_C-plus-plus_Reference/) is used to optimize the hRNG
- to use the optimization feature, specificy the install prefix within the main CMakeLists.txt file
- install prefix as the SET(NLOPT_HOME ".../nlopt/")" to where nlopt is installed locally
