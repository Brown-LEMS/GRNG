#include "optimization.hpp"

#include <nlopt.hpp>

#include "hrng-data.hpp"
#include "incremental-build.hpp"
#include "online-search.hpp"
#include "results-saver.hpp"

Optimization::Optimization(unsigned int const dimension, float *const &dataPointer, unsigned int const datasetSize, float *const &testPointer,
                           unsigned int const testsetSize, int const numThreads, std::string luneType) {
    _data = hRNGData(dimension, dataPointer, datasetSize, testPointer, testsetSize, numThreads, luneType);
}

Optimization::Optimization(unsigned int const dimension, float *const &dataPointer, unsigned int const datasetSize, float *const &testPointer,
                           unsigned int const testsetSize, int const numThreads, std::string luneType, bool cacheAll) {
    _data = hRNGData(dimension, dataPointer, datasetSize, testPointer, testsetSize, numThreads, luneType, cacheAll);
}


/**
 * @brief Functions to help save opt. info to output file
 *
 */
namespace ResultsSaver_Optimization {
void printResultsHeader(ResultsSaver &saver);
void printResults(ResultsSaver &saver, hRNGResults &results);
void printBestResults(ResultsSaver &saver, hRNGResults &results);
}  // namespace ResultsSaver_Optimization

/**
 * @brief Objective function to minimize for NLOPT. runs hRNG, looking to optimize search distance computations
 *
 * @param x
 * @param grad
 * @param hRNGData_Void
 * @return double
 */
double Optimization::objectiveFunction(const std::vector<double> &x, std::vector<double> &grad, void *hRNGData_Void) {
     // convert void pointer back to hRNG_Data
    hRNGData *data = (hRNGData *)hRNGData_Void;
    data->_results._buildCount++;

    // create effective pivot radius vector from double to float
    std::vector<float> pivotRadiusVector(x.begin(), x.end());
    pivotRadiusVector.push_back(0.0f);
    data->_results._pivotRadiusVector = pivotRadiusVector;

    // effective pivot radius vector, convert to float from double
    printf("hRNG(");
    std::vector<float> effectivePivotRadiusVector = pivotRadiusVector;
    for (int i = 0; i < (int)effectivePivotRadiusVector.size(); i++) {
        for (int j = i + 1; j < (int)x.size(); j++) {
            effectivePivotRadiusVector[i] += (float)pivotRadiusVector[j];
        }
        printf("%.6f,", pivotRadiusVector[i]);  // print radii value
    }
    printf(") = ");
    data->_results._effectivePivotRadiusVector = effectivePivotRadiusVector;


    // pivot layers
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectivePivotRadiusVector, data->_luneType, data->_dataPointer, data->_datasetSize, data->_dimension,data->_cacheAll);

    // incremnetal build
    IncrementalBuild IB(pivotLayers);
    IB.setNumThreads(data->_numThreads);
    IB.build(true);

    data->_results._numberOfPivotsVector = IB._numberOfPivotsPerLayer;
    data->_results._buildDistances = IB._distanceCount_build;
    data->_results._buildTime = IB._time_build.count();
    data->_results._averageDistances = IB._distanceCount_averageQuery;
    data->_results._averageTime = IB._time_averageQuery.count();
    data->_results._memoryUsage = (float)(IB._currentRSS / (double)(1 << 30));
    data->_results._linkCountVector = IB._numberLinksPerLayer;
    data->_results._averageDegreeVector = IB._averageDegreePerLayer;
    data->_results._minDegreeVector = IB._minDegreePerLayer;
    data->_results._maxDegreeVector = IB._maxDegreePerLayer;

    // online search
    OnlineSearch OS(pivotLayers);
    OS.setNumThreads(data->_numThreads);
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> RNG_neighors{};
    OS.search(data->_testPointer, data->_testsetSize, data->_dimension, true, RNG_neighors);

    data->_results._searchDistances = OS._distanceCount_averageQuery;
    data->_results._searchTime = OS._time_averageQuery.count() * 1000;  // to milliseconds

    // save results to results saver
    // data->resultsSaver_saveResults();
    // data->resultsSaver_keepBest();
    ResultsSaver_Optimization::printResults(*data->_saver, data->_results);

    printf("Search: %.3f, Incremental: %.3f\n", data->_results._searchDistances, data->_results._averageDistances);
    return (double) data->_results._searchDistances; //data->_results._averageDistances; //
}

bool Optimization::localOptimization(std::string method, unsigned n, std::vector<double> &x, double &f, std::vector<double> lowerBounds,
                                     std::vector<double> upperBounds, std::vector<double> stepSize, std::vector<double> xTolVector, double fTol,
                                     int maxEval, std::string resultsSaverPath) {
    printf("Using NLOPT Local Optimization Functions...\n");
    bool success = true;

    // prepare results saver
    _data._saver = new ResultsSaver(resultsSaverPath);
    ResultsSaver_Optimization::printResultsHeader(*_data._saver);

    // create description of optimization parameters
    // std::string test_description = std::string("Optimization Test: ")
    //                                    .append("; D; ")
    //                                    .append(std::to_string(data_.dimension_))
    //                                    .append("; N; ")
    //                                    .append(std::to_string(data_.datasetSize_))
    //                                    .append("; lune; ")
    //                                    .append(data_.luneType_)
    //                                    .append("; L; ")
    //                                    .append(std::to_string(n + 1))
    //                                    .append("; Local Algo; ")
    //                                    .append(method);

    // strings for each vector to append to description
    // std::string sr, lb, ub, is, xT;
    // for (int i = 0; i < (int)n; i++) {
    //     sr = sr.append(std::to_string(x[i])).append(" ");
    //     lb = lb.append(std::to_string(lowerBounds[i])).append(" ");
    //     ub = ub.append(std::to_string(upperBounds[i])).append(" ");
    //     is = lb.append(std::to_string(stepSize[i])).append(" ");
    //     if (xTolVector.size() > 0) xT = xT.append(std::to_string(xTolVector[i])).append(" ");
    // }
    // test_description = test_description.append("; startingRadii; ")
    //                        .append(sr)
    //                        .append("; lowerBounds; ")
    //                        .append(lb)
    //                        .append("; upperBounds; ")
    //                        .append(ub)
    //                        .append("; fTol; ")
    //                        .append(std::to_string(fTol))
    //                        .append("; stepSize; ")
    //                        .append(is)
    //                        .append("; xTolVector; ")
    //                        .append(xT)
    //                        .append("; maxEval; ")
    //                        .append(std::to_string(maxEval));

    // save test description and header for results
    // std::cout << test_description << std::endl;
    // saver.printLine(test_description);
    // data_.resultsSaver_saveHeader(n + 1);

    // convert hRNG_Data to void* for NLOPT
    // hRNGData _data;
    void *hRNG_data = (void *)&_data;

    // create optimization object
    // add bounds as algorithms allow
    nlopt::opt opt;
    if (method == "LN_COBYLA") {
        // COBYLA (Constrained Optimization BY Linear Approximations)
        opt = nlopt::opt(nlopt::LN_COBYLA, n);
        opt.set_lower_bounds(lowerBounds);  // sets the lower bounds
        opt.set_upper_bounds(upperBounds);  // sets the upper bounds
    } else {
        printf("Unrecognized Optimization Method...\n");
        return 0;
    }

    // add initial stepSize
    opt.set_initial_step(stepSize);

    // set objective (Minimize hRNG), passing our information
    opt.set_min_objective(objectiveFunction, hRNG_data);

    // stopping criteria via relative change in objective value
    opt.set_ftol_rel(fTol);

    // stopping criteria via additive step size
    if (xTolVector.size() > 0) {
        opt.set_xtol_abs(xTolVector);
    }

    // max evalutations
    opt.set_maxeval(maxEval);

    // perform optimization
    printf("Performing Optimization...\n");
    std::chrono::high_resolution_clock::time_point tStart = std::chrono::high_resolution_clock::now();
    int resultInt = opt.optimize(x, f);
    std::chrono::high_resolution_clock::time_point tEnd = std::chrono::high_resolution_clock::now();
    float optimizationTime = (float)(std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart)).count();
    //_data.resultsSaver_saveBestResults(optimizationTime);
    _data._saver->saveVectorElement(x);
    _data._saver->saveElement(f);
    delete _data._saver;

    // print/save optimization termination reasoning
    std::string resultDescription = std::string("Terminated by Command ").append(std::to_string(resultInt));
    switch (resultInt) {
    case 1:
        // saver.printLine(resultDescription.append(": Success."));
        printf("%s: success!\n", resultDescription.c_str());
        break;
    case 2:
        // saver.printLine(resultDescription.append(": stopVal Reached."));
        printf("%s: stopVal Reached!\n", resultDescription.c_str());
        break;
    case 3:
        // saver.printLine(resultDescription.append(": fTol Reached."));
        printf("%s: fTol Reached!\n", resultDescription.c_str());
        break;
    case 4:
        // saver.printLine(resultDescription.append(": xTol Reached."));
        printf("%s: xTol Reached!\n", resultDescription.c_str());
        break;
    case 5:
        // saver.printLine(resultDescription.append(": maxEval Reached."));
        printf("%s: maxEval Reached!\n", resultDescription.c_str());
        break;
    case 6:
        // saver.printLine(resultDescription.append(": maxTime Reached."));
        printf("%s:  maxTime Reached!\n", resultDescription.c_str());
        break;
    case -1:
        // saver.printLine(resultDescription.append(": Generic Code Failure."));
        printf("%s: Generic Code Failure!\n", resultDescription.c_str());
        success = false;
        break;
    case -2:
        // saver.printLine(resultDescription.append(": Invalid Arguments."));
        printf("%s: Invalid Arguments!\n", resultDescription.c_str());
        success = false;
        break;
    case -3:
        // saver.printLine(resultDescription.append(": Ran out of Memory."));
        printf("%s: Ran out of Memory!\n", resultDescription.c_str());
        success = false;
        break;
    case -4:
        // saver.printLine(resultDescription.append(": Roundoff Errors."));
        printf("%s: Roundoff Errors!\n", resultDescription.c_str());
        success = false;
        break;
    case -5:
        // saver.printLine(resultDescription.append(": Forced Termination."));
        printf("%s: Forced Termination!\n", resultDescription.c_str());
        success = false;
        break;
    default:
        break;
    }
    return success;
}

void ResultsSaver_Optimization::printResultsHeader(ResultsSaver &saver) {
    std::string del(1, saver._del);
    std::string header = "";
    header = header.append("Build").append(del);
    header = header.append("Pivot Radius Vector").append(del);
    header = header.append("Effective Radii Vector").append(del);
    header = header.append("Number of Pivots Vector").append(del);
    header = header.append("Build Distances").append(del);
    header = header.append("Build Time (s)").append(del);
    header = header.append("Average Search Distances").append(del);
    header = header.append("Average Search Time (ms)").append(del);
    header = header.append("Average Incremental Distances").append(del);
    header = header.append("Average Incremental Time (ms)").append(del);
    header = header.append("PeakRSS (GB)").append(del);
    header = header.append("Link Count Vector").append(del);
    header = header.append("Average Degree Vector").append(del);
    header = header.append("Min Degree Vector").append(del);
    header = header.append("Max Degree Vector");
    saver.printLine(header);
    return;
}

void ResultsSaver_Optimization::printResults(ResultsSaver &saver, hRNGResults &results) {
    std::string del(1, saver._del);
    std::string header = "";
    saver.saveElement(results._buildCount);
    saver.saveVectorElement(results._pivotRadiusVector);
    saver.saveVectorElement(results._effectivePivotRadiusVector);
    saver.saveVectorElement(results._numberOfPivotsVector);
    saver.saveElement(results._buildDistances);
    saver.saveElement(results._buildTime);
    saver.saveElement(results._searchDistances);
    saver.saveElement(results._searchTime);
    saver.saveElement(results._averageDistances);
    saver.saveElement(results._averageTime);
    saver.saveElement(results._memoryUsage);
    saver.saveVectorElement(results._linkCountVector);
    saver.saveVectorElement(results._averageDegreeVector);
    saver.saveVectorElement(results._minDegreeVector);
    saver.saveVectorElement(results._maxDegreeVector);
    saver.newLine();
    return;
}

void printResults(ResultsSaver &saver);