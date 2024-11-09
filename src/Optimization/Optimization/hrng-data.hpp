#ifndef hRNGData_hpp
#define hRNGData_hpp

#include "hrng-results.hpp"
#include "results-saver.hpp"

/**
 * @brief Struct to hold information during optimization. Passed as void* to nlopt. keeps track of hRNG parameters and holds results. Also saves
 * results
 *
 */
struct hRNGData {

    // constructor
    hRNGData(){};
    hRNGData(unsigned int const dimension, float* const& dataPointer, unsigned int const datasetSize, float* const& testPointer,
             unsigned int const testsetSize, int numThreads, std::string luneType)
        : _dimension(dimension),
          _dataPointer(dataPointer),
          _datasetSize(datasetSize),
          _testPointer(testPointer),
          _testsetSize(testsetSize),
          _numThreads(numThreads),
          _luneType(luneType) {};
    
    hRNGData(unsigned int const dimension, float* const& dataPointer, unsigned int const datasetSize, float* const& testPointer,
             unsigned int const testsetSize, int numThreads, std::string luneType, bool cacheAll)
        : _dimension(dimension),
          _dataPointer(dataPointer),
          _datasetSize(datasetSize),
          _testPointer(testPointer),
          _testsetSize(testsetSize),
          _numThreads(numThreads),
          _luneType(luneType),
          _cacheAll(cacheAll) {};

    ~hRNGData(){};

    // data parameters
    unsigned int _buildCount = 0;
    int _numLayers = 0;

    // results holders
    hRNGResults _results;
    hRNGResults _bestResults;
    ResultsSaver* _saver;

    // build and search parameters
    unsigned int _dimension = 0;
    float* _dataPointer = NULL;
    unsigned int _datasetSize = 0;
    float* _testPointer = NULL;
    unsigned int _testsetSize = 0;
    int _numThreads = 0;
    std::string _luneType = "";
    bool _cacheAll = false;
};

#endif  // hRNGData_hpp