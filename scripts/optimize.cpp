#include <chrono>
#include <cstdio>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "optimization.hpp"
#include <string>

std::string getDateTime() {
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d", std::localtime(&now));
    return s;
}

int main(int argc, char** argv) {
    printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    unsigned int datasetSize = 1600;
    unsigned int dataSeed = 3;
    unsigned int testsetSize = 100;
    unsigned int testSeed = 5;
    std::string luneType = "Berk";
    int numThreads = 1;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size Vector");
    app.add_option("-T,--testsetSize", testsetSize, "Dataset Size Vector");
    app.add_option("-l,--luneType", luneType, "Lune Type [Berk, Cole]");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");

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

    float* dataPointer = NULL;
    float* testPointer = NULL;
    if (dataset == "uniform") {
        Datasets::uniformDataset(dataPointer, dimension, datasetSize, dataSeed);
        Datasets::uniformDataset(testPointer, dimension, testsetSize, testSeed);
    } else {
        printf("Unrecognized dataset!\n");
        return 0;
    }

    // initialize optimizer
    Optimization optimize(dimension, dataPointer, datasetSize, testPointer, testsetSize, numThreads, luneType);

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
    std::string outputDirectory = std::string("results/") + getDateTime().c_str() + "/";
    printf("Created Folder For Results: %s\n", outputDirectory.c_str());
    if (mkdir(outputDirectory.c_str(), ACCESSPERMS) == 0) {
        printf("Created Folder For Results: %s\n", outputDirectory.c_str());
    }
    std::string output_path = std::string(dataset)
                                  .append("_D-")
                                  .append(std::to_string(dimension))
                                  .append("_N-")
                                  .append(std::to_string(datasetSize))
                                  .append("_L-")
                                  .append(std::to_string(x.size() + 1))
                                  .append("_")
                                  .append(luneType)
                                  .append("_n-")
                                  .append(std::to_string(numThreads)).append(".txt");

    std::string out = outputDirectory.append(output_path);
    printf("Output Path: %s\n", out.c_str());

    // perform optimization
    optimize.localOptimization(method, n, x, f, lowerBoundsVector, upperBoundsVector, stepSizeVector, xTolVector, fTol, maxEval, out);

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    delete[] testPointer;
    return 0;
}
