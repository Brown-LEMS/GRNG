#define BOOST_TEST_DYN_LINK  // deals with no main() issue
#define BOOST_TEST_MODULE OptimizationTest
#include <boost/test/unit_test.hpp>

#include "datasets.hpp"
#include "optimization.hpp"

BOOST_AUTO_TEST_SUITE(optimization)

BOOST_AUTO_TEST_CASE(neighbors) {
    unsigned int const dimension = 2;

    // create dataset
    unsigned int const datasetSize = 1600;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);

    // create testset
    unsigned int const testsetSize = 100;
    float* testPointer;
    Datasets::uniformDataset(testPointer, dimension, testsetSize, 3);

    std::string method = "LN_COBYLA";
    unsigned int n = 3;
    std::vector<double> x = {0.21872, 0.08995, 0.02123};
    double f = 0;
    std::vector<double> lowerBounds = {0, 0, 0};
    std::vector<double> upperBounds = {100, 100, 100};
    std::vector<double> stepSize = {0.1, 0.1, 0.1};
    std::vector<double> xTolVector = std::vector<double>{};
    double fTol = 0.1;
    int maxEval = 20;
    std::string resultsSaverPath = "./here.txt";

    // make a class
    Optimization optimizer(dimension, dataPointer, datasetSize, testPointer, testsetSize, 1, "Berk");
    optimizer.localOptimization(method, n, x, f, lowerBounds, upperBounds, stepSize, xTolVector, fTol, maxEval, resultsSaverPath);

    delete[] dataPointer;
    delete[] testPointer;
    return;
}

BOOST_AUTO_TEST_SUITE_END()