#include <chrono>
#include <cstdio>
#include <string>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "optimization.hpp"

std::string getDateTime() {
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d", std::localtime(&now));
    return s;
}

int main(int argc, char** argv) {
    printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    unsigned int startSize = 400;
    unsigned int endSize = 1600;
    unsigned int dataSeed = 3;
    unsigned int testsetSize = 100;
    unsigned int testSeed = 5;
    std::string luneType = "Berk";
    int numThreads = 1;
    bool cacheAll = false;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-L,--startSize", startSize, "Dataset Size");
    app.add_option("-N,--endSize", endSize, "Dataset Size");
    app.add_option("-T,--testsetSize", testsetSize, "testset Size");
    app.add_option("-l,--luneType", luneType, "Lune Type [Berk, Cole]");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-C,--cacheAll", cacheAll, "Cache all distances");

    // optimization parameters
    std::vector<float> startingRadii{};
    std::vector<float> stepSize{};
    std::vector<float> lowerBounds{};
    std::vector<float> upperBounds{};
    float fTol = 0.01;
    int maxEval = 1000;
    app.add_option("-r,--startingRadii", startingRadii, "radiusVector for hierarchy (leave out RNG)")->required();
    app.add_option("-s,--stepSize", stepSize, "radiusVector for hierarchy (leave out RNG)");
    app.add_option("-[,--lowerBounds", lowerBounds, "upperBounds of radii");
    app.add_option("-],--upperBounds", upperBounds, "lowerBounds of radii");
    app.add_option("-f,--fTol", fTol, "Stopping Criteria: tolerance for change in f");
    app.add_option("-M,--maxEval", maxEval, "Stopping Criteria: max function evaluations");

    CLI11_PARSE(app, argc, argv);
    unsigned int numLayers = startingRadii.size();
    unsigned int datasetSize = endSize;

    // get full dataset, in terms of endsize
    float* dataPointer = NULL;
    float* testPointer = NULL;
    float* validationPointer = NULL;
    if (dataset == "uniform") {
        // if uniform, getting full endSize then getting testset
        Datasets::uniformDataset(dataPointer, dimension, datasetSize, dataSeed);
        Datasets::uniformDataset(validationPointer, dimension, testsetSize, testSeed);
        Datasets::uniformDataset(testPointer, dimension, testsetSize, testSeed);

    } else if (dataset == "cluster") {
        Datasets::clusterDataset(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer, dimension, datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize);  // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "Cities") {
        Datasets::Cities(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer, dimension, datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize);  // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "Corel1") {
        Datasets::Corel68k1(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "Corel2") {
        Datasets::Corel68k2(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "Corel3") {
        Datasets::Corel68k3(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "SIFT") {
        std::string data_path = std::string(dataDirectory).append("sift/sift_base.fvecs");
        Datasets::SIFT1M(dataPointer, dimension, datasetSize, data_path);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "MNIST") {
        Datasets::MNIST(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "DEEP") {
        Datasets::DEEP10M(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // reserved
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation
    } else if (dataset == "LA") {
        Datasets::LA(dataPointer, dimension, datasetSize,testPointer,testsetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer, dimension, datasetSize);
        Datasets::extractTestset(dimension, dataPointer, datasetSize, validationPointer, testsetSize);  // validation

    } else {
        printf("Unrecognized dataset!\n");
        return 0;
    }

    // initialize optimization parameters
    std::string method = "LN_COBYLA";
    unsigned int n = startingRadii.size();
    std::vector<double> x(startingRadii.begin(), startingRadii.end());  // returns the best radii
    double f = 0;                                                       // returns the best optimization value
    std::vector<double> lowerBoundsVector(n, 0);
    if (lowerBounds.size() > 0) {
        lowerBoundsVector = std::vector<double>(lowerBounds.begin(), lowerBounds.end());
    }
    std::vector<double> upperBoundsVector(n, 10000);
    if (upperBounds.size() > 0) {
        upperBoundsVector = std::vector<double>(upperBounds.begin(), upperBounds.end());
    }
    std::vector<double> stepSizeVector(n, 0.01);
    if (stepSize.size() > 0) {
        stepSizeVector = std::vector<double>(stepSize.begin(), stepSize.end());
    }
    std::vector<double> xTolVector(n, 0.00001);

    // create output path
    //std::string resultsHomeDirectory = "/home/colef8/Documents/research/hRNG/SISAP/bin/results/";
    std::string resultsHomeDirectory = "/users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/results/";
    std::string resultsDateDirectory = std::string(resultsHomeDirectory).append(getDateTime().c_str()).append("/");
    mkdir(resultsDateDirectory.c_str(), ACCESSPERMS);
    std::string resultsOptDirectory = std::string(resultsDateDirectory).append(dataset).append("/");
    mkdir(resultsOptDirectory.c_str(), ACCESSPERMS);
    std::string resultsDirectory = std::string(resultsOptDirectory).append("optimization/");
    mkdir(resultsDirectory.c_str(), ACCESSPERMS);

    // begin opt loop
    unsigned int optCount = 0;
    unsigned int optSize = startSize;
    while (optSize <= endSize) {
        optCount++;

        // initialize optimizer
        Optimization optimize(dimension, dataPointer, optSize, validationPointer, testsetSize, numThreads, luneType, cacheAll);

        std::string output_path = std::string(dataset)
                                      .append("_D-")
                                      .append(std::to_string(dimension))
                                      .append("_L-")
                                      .append(std::to_string(numLayers + 1))
                                      .append("_N-")
                                      .append(std::to_string(optSize))
                                      .append("_OC-")
                                      .append(std::to_string(optCount))
                                      .append("_")
                                      .append(luneType)
                                      .append("_n-")
                                      .append(std::to_string(numThreads))
                                      .append(".txt");

        std::string out = std::string(resultsDirectory).append(output_path);
        printf("Output Path: %s\n", out.c_str());

        // perform optimization
        optimize.localOptimization(method, n, x, f, lowerBoundsVector, upperBoundsVector, stepSizeVector, xTolVector, fTol, maxEval, out);

        // update parameters for next round
        optSize *= 2;
        for (int i = 0; i < x.size(); i++) {
            x[i] -= stepSizeVector[i];
        }
    }

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    delete[] testPointer;
    delete[] validationPointer;
    return 0;
}