#include <cstdio>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "offline-build.cpp"

std::vector<float> getEffectiveRadiusVector(std::vector<float> const vec)
{
    std::vector<float> new_vec = vec;
    for (int i = 0; i < (int)vec.size(); i++)
    {
        for (int j = i + 1; j < (int)vec.size(); j++)
        {
            new_vec[i] += (float)vec[j];
        }
    }
    return new_vec;
}

int main(int argc, char **argv)
{
    printf("Begin Index Construction... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 0;
    unsigned int datasetSize = 0;
    unsigned int dataSeed = 3;
    unsigned int testsetSize = 100;
    unsigned int testSeed = 5;
    std::string luneType = "Berk";
    std::vector<float> radiusVector{}; // REQUIRED
    bool effectiveRadius = true;
    bool cacheAll = false;
    int numThreads = 1;
    bool verbose = true;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size Vector");
    app.add_option("-T,--testsetSize", testsetSize, "Dataset Size Vector");
    app.add_option("-l,--luneType", luneType, "Lune Type [Berk, Cole]");
    app.add_option("-C,--cacheAll", cacheAll, "cache all distance computations, or just relationships");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-r,--radiusVector", radiusVector, "radiusVector for hierarchy (leave out RNG)")->required();
    app.add_option("-e,--effectiveRadius", effectiveRadius, "Use effective radius size");
    app.add_option("-v,--verbose", verbose, "Show hRNG Statistics During Eval");
    CLI11_PARSE(app, argc, argv);

    // get dataset
    float *dataPointer = NULL;
    float *testPointer = NULL;
    if (dataset == "uniform")
    {
        Datasets::uniformDataset(dataPointer, dimension, datasetSize, dataSeed);
        Datasets::uniformDataset(testPointer, dimension, testsetSize, testSeed);
    }
    else
    {
        printf("Unrecognized dataset: %s\n", dataset.c_str());
        return 0;
    }
    printf("Dataset: %s, Dimension: %u, DatasetSize: %u\n", dataset.c_str(), dimension, datasetSize);
    printf("TestsetSize: %u\n", testsetSize);

    // get effective radius vector
    radiusVector.push_back(0);
    std::vector<float> effectivePivotRadiusVector = getEffectiveRadiusVector(radiusVector);

    // initialize hierarchy
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectivePivotRadiusVector, luneType, dataPointer, datasetSize, dimension, true);

    // incrementally build hierarchy
    OfflineBuild OB(pivotLayers);
    OB.selectPivots();

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    delete[] testPointer;
    return 0;
}
