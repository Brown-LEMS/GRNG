#include <cstdio>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "hierarchical-rng.hpp"

int main(int argc, char** argv) {
    printf("Begin Index Construction... \n");
    CLI::App app{"Index Construction"};
    unsigned int dimension;    // REQUIRED
    unsigned int datasetSize;  // REQUIRED
    unsigned int dataSeed = 3;
    std::string luneType = "Berk";
    std::vector<float> radiusVector;  // REQUIRED
    bool cacheAll = true;
    bool incremental = false;
    int numThreads = 1;
    bool verbose = false;
    app.add_option("-D,--dimension", dimension, "Dimension")->required();
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size Vector")->required();
    app.add_option("-l,--luneType", luneType, "Lune Type [Berk, Cole]");
    app.add_option("-C,--cacheAll", cacheAll, "cache all distance computations, or just relationships");
    app.add_option("-i,--incremental", incremental, "Incremental or Batch build");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-r,--radiusVector", radiusVector, "radiusVector for hierarchy (leave out RNG)")->required();
    app.add_option("-v,--verbose", verbose, "Show hRNG Statistics During Eval");
    CLI11_PARSE(app, argc, argv);

    // get dataset
    float* dataPointer = NULL;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, dataSeed);

    // initialize hRNG
    HierarchicalRNG hRNG(dataPointer, datasetSize, dimension, cacheAll);
    radiusVector.push_back(0);
    std::vector<float> effectivePivotRadiusVector = hRNG.getEffectiveRadiusVector(radiusVector);
    if (incremental) {
        hRNG.incrementalBuild(effectivePivotRadiusVector, luneType, numThreads, verbose);
    } else {
        hRNG.offlineBuild(effectivePivotRadiusVector, luneType, numThreads, verbose);
    }

    // should have an option to verify correct RNG by brute force comparison

    printf("Done! Have a good day! \n");
    return 0;
}