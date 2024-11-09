// Cole Foster

#include <chrono>
#include <cstdio>
#include <vector>

#include "CLI11.hpp"
#include "datasets.hpp"
// #include "hierarchical-rng.hpp"

#include "incremental-build.hpp"
#include "online-search.hpp"
#include "pivot-layer-hierarchy.hpp"
#include "verification.hpp"
#include "memory-usage.h"
#include "rayar.hpp"

class ResultsOuputter {
  public:
    ResultsOuputter(){};
    ResultsOuputter(std::string filename) : _filename(filename) {
        printf("Initiated File Saver for: %s\n", _filename.c_str());
        _getDateTime();  // retrieve datetime

        // looping through filenames
        if (_fileExists(_filename)) {
            printf("Whoops! filename: '%s' exists, appending date and time to end...\n", _filename.c_str());
            _filename = std::string(filename).append(_dateTime);
            printf("New filename: %s\n", _filename.c_str());
        }

        // initialize
        _file.open(_filename);
        _file.close();
    };

    ~ResultsOuputter(){};

    // member vars
    std::string _filename = "";
    std::string _dateTime = "";
    std::ofstream _file;
    char _del = ',';
    char _vecDel = ';';

    inline void _getDateTime() {
        std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string s(30, '\0');
        std::strftime(&s[0], s.size(), "_%Y-%m-%d_%H:%M:%S", std::localtime(&now));
        _dateTime = s;
        return;
    }

    inline bool _fileExists(const std::string& filename) {
        struct stat buf;
        if (stat(filename.c_str(), &buf) != -1) {
            return true;
        }
        return false;
    }

    inline void printLine(std::string const line) {
        _file.open(_filename.c_str(), std::ios_base::app);  // append to file, instead of overwrite
        _file << line.c_str() << '\n';
        _file.close();
    }

    inline void printDateTime() { printLine(_dateTime); }

    // saves vector, as an entire new line, '\n' after
    template <class T>
    inline void saveVectorLine(std::vector<T> vec) {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        for (int i = 0; i < (int)vec.size(); i++) _file << vec[i] << _del;
        _file << '\n';
        _file.close();
    }

    // saves vector using delimiter of ';' and wrapped by { }
    // ASSUMES FILE IS ALREADY OPEN
    template <class T>
    inline void saveVectorElement(std::vector<T> vec) {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        _file << "{" << _vecDel;
        for (int i = 0; i < (int)vec.size(); i++) _file << vec[i] << _vecDel;
        _file << "}" << _del;
        _file.close();
    }

    template <class T>
    inline void saveElement(T element) {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        _file << element << _del;
        _file.close();
    }

    inline void newLine() {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        _file << '\n';
        _file.close();
    }
};

std::string getDateTime() {
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d", std::localtime(&now));
    return s;
}

std::vector<float> getEffectiveRadiusVector(std::vector<float> const vec) {
    std::vector<float> new_vec = vec;
    for (int i = 0; i < (int)vec.size(); i++) {
        for (int j = i + 1; j < (int)vec.size(); j++) {
            new_vec[i] += (float)vec[j];
        }
    }
    return new_vec;
}

int main(int argc, char** argv) {
    printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 0;
    unsigned int datasetSize = 0;
    unsigned int dataSeed = 3;
    unsigned int testsetSize = 100;
    unsigned int testSeed = 5;
    int numThreads = 1;
    float epsilon = 0.1;
    int L = 4;
    bool GT = false;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size Vector");
    app.add_option("-T,--testsetSize", testsetSize, "Dataset Size Vector");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-L,--Lorder", L, "Rayar Neighbor Order");
    app.add_option("-G,--GT", GT, "Naive RNG comparison?");
    CLI11_PARSE(app, argc, argv);

    float* dataPointer = NULL;
    float* testPointer = NULL;
    if (dataset == "uniform") {
        Datasets::uniformDataset(dataPointer, dimension, datasetSize, dataSeed);
        Datasets::uniformDataset(testPointer, dimension, testsetSize, testSeed);
    } else if (dataset == "Cities") {
        Datasets::Cities(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer, dimension, datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize);  // queries
    } else if (dataset == "Corel1") {
        Datasets::Corel68k1(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // queries
    } else if (dataset == "Corel2") {
        Datasets::Corel68k2(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // queries
    } else if (dataset == "Corel3") {
        Datasets::Corel68k3(dataPointer, dimension, datasetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // queries
    } else if (dataset == "SIFT") {
        std::string data_path = std::string(dataDirectory).append("sift/sift_base.fvecs");
        Datasets::SIFT1M(dataPointer, dimension, datasetSize, data_path);
        Datasets::datasetShuffle(dataPointer,dimension,datasetSize);
        testsetSize = 100;
        Datasets::extractTestset(dimension, dataPointer, datasetSize, testPointer, testsetSize); // queries
    } else if (dataset == "LA") {
        Datasets::LA(dataPointer, dimension, datasetSize,testPointer,testsetSize, dataDirectory);
        Datasets::datasetShuffle(dataPointer, dimension, datasetSize);

    } else {
        printf("Unrecognized dataset: %s\n", dataset.c_str());
        return 0;
    }

    std::string resultsHomeDirectory = "/users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/results/";
    std::string resultsDateDirectory = std::string(resultsHomeDirectory).append(getDateTime().c_str()).append("/");
    mkdir(resultsDateDirectory.c_str(), ACCESSPERMS);
    std::string resultsDirectory = std::string(resultsDateDirectory).append(dataset).append("/");
    mkdir(resultsDirectory.c_str(), ACCESSPERMS);
    std::string resultsName = std::string(dataset)
                                  .append("_D-")
                                  .append(std::to_string(dimension))
                                  .append("_N-")
                                  .append(std::to_string(datasetSize))
                                  .append("_RAYAR_L-4.txt");
    std::string resultsPath = std::string(resultsDirectory).append(resultsName);
    printf("Saving Results to: %s \n", resultsPath.c_str());

    ResultsOuputter results(resultsPath);
    results.printDateTime();
    results.saveElement(dimension);
    results.saveElement(datasetSize);
    results.saveElement(L);
    results.saveElement(numThreads);
    results.newLine();


    // perform index construction
    std::vector<tsl::sparse_set<unsigned int>> RNG_neighbors4{};
    unsigned long long int build_distanceComputations4 = 0;
    float build_time4 = 0;
    rayar::approximateRNG(dataPointer, datasetSize, dimension, epsilon, L, numThreads, RNG_neighbors4, build_distanceComputations4, build_time4);

    results.saveElement(build_distanceComputations4);
    results.saveElement(build_time4);
    results.saveElement((double)((double) getCurrentRSS() / (double)(1 << 30)));

    // perform GT
    std::vector<tsl::sparse_set<unsigned int>> GT_RNG_neighbors{};
    unsigned long long int GT_build_distanceComputations = 0;
    float GT_build_time = 0;
    if (GT) {
        printf("Using Naive RNG Construction...\n");

        // get effective radius vector
        std::vector<float> radiusVector = {489.518,122.797,30.8038,7.72721};
        radiusVector.push_back(0);
        std::vector<float> effectivePivotRadiusVector = getEffectiveRadiusVector(radiusVector);

        // initialize hierarchy
        std::shared_ptr<PivotLayerHierarchy> pivotLayers =
            std::make_shared<PivotLayerHierarchy>(effectivePivotRadiusVector, "Berk", dataPointer, datasetSize, dimension, false);

        // incrementally build hierarchy
        IncrementalBuild IB(pivotLayers);
        IB.setNumThreads(numThreads);
        IB.build(true);

        tsl::sparse_map<unsigned int,tsl::sparse_set<unsigned int>> neighbors1{};
        pivotLayers->getRNG(neighbors1);

        GT_RNG_neighbors.resize(neighbors1.size());
        tsl::sparse_map<unsigned int,tsl::sparse_set<unsigned int>>::const_iterator it1;
        for (it1 = neighbors1.begin(); it1 != neighbors1.end(); it1++) {
            unsigned int const index = (*it1).first;
            tsl::sparse_set<unsigned int> const& neighbs = (*it1).second;
            GT_RNG_neighbors[index].insert(neighbs.begin(),neighbs.end());
        }
    }
    results.saveElement(GT_build_distanceComputations);
    results.saveElement(GT_build_time);

    // perform search by rayar
    std::vector<tsl::sparse_set<unsigned int>> RNG_searchNeighbors4{};
    double search_distanceComputations4 = 0;
    float search_time4 = 0;
    rayar::search(dimension, dataPointer, datasetSize, testPointer, testsetSize, epsilon, numThreads, RNG_neighbors4,
                  RNG_searchNeighbors4, search_distanceComputations4, search_time4);
    results.saveElement(search_distanceComputations4);
    results.saveElement(search_time4);
    

    // compute brute force neighbors
    std::vector<tsl::sparse_set<unsigned int>> GT_RNG_searchNeighbors{};
    float GT_search_distanceComputations = 0;
    float GT_search_time = 0;
    {                 
        verification::bruteForceSearch(dimension, dataPointer, datasetSize, testPointer, testsetSize, GT_RNG_searchNeighbors,
                                       GT_search_distanceComputations, GT_search_time);
    }
    results.saveElement(GT_search_distanceComputations);
    results.saveElement(GT_search_time);

    
    // accurarcy of main rng 
    bool correct;
    unsigned int numberOfLinks = 0, extraLinks = 0, missingLinks = 0;
    if (GT) {
        correct = verification::graphComparison(GT_RNG_neighbors, RNG_neighbors4, numberOfLinks, extraLinks, missingLinks);
        printf("Rayar4 RNG | Distances, Time, Correct, Extra, Missing: \n%u\n%.3f\n%u\n%u\n%u\n", build_distanceComputations4, build_time4, correct, extraLinks, missingLinks);
    }
    results.saveElement(correct);
    results.saveElement(numberOfLinks);
    results.saveElement(extraLinks);
    results.saveElement(missingLinks);

    correct = verification::neighborsComparison(GT_RNG_searchNeighbors, RNG_searchNeighbors4, numberOfLinks, extraLinks, missingLinks);
    results.saveElement(correct);
    results.saveElement(numberOfLinks);
    results.saveElement(extraLinks);
    results.saveElement(missingLinks);

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    delete[] testPointer;
    return 0;
}
