#include "Test_IncrementalBuild.hpp"
#include "offline-build.hpp"
#include "datasets.hpp"

std::vector<float> getEffectiveRadiusVector(std::vector<float> vec) {
    std::vector<float> new_vec = vec;
    for (int i = 0; i < (int)vec.size(); i++) {
        for (int j = i + 1; j < (int)vec.size(); j++) {
            new_vec[i] += (float)vec[j];
        }
    }
    return new_vec;
}

void printMapSet(tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> mapSet) {
    std::cout << "{ \n";
    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>::const_iterator it1;
    for (it1 = mapSet.begin(); it1 != mapSet.end(); it1++) {
        unsigned int key = (*it1).first;
        std::cout << "  " << key << ":{";
        tsl::sparse_set<unsigned int> set = (*it1).second;
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = set.begin(); it1 != set.end(); it1++) {
            std::cout << (*it1) << ",";
        }
        std::cout << "}, \n";
    }
    std::cout << "}" << std::endl;
};

void printSet9(tsl::sparse_set<unsigned int> set) {
    std::cout << "{";
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        std::cout << (*it1) << ",";
    }
    std::cout << "}, \n";
};



void Test_IncrementalBuild::test_pivotSelection() {
    std::cout << "Unit Test: Pivot Selection" << std::endl;

    // create dataset
    unsigned int const datasetSize = 800;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0.11395, 0.04995, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers1 =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension);

    std::shared_ptr<PivotLayerHierarchy> pivotLayers2 =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension);

    // offline build to get RNG
    OfflineBuild OB(pivotLayers1);
    OB.build(false);

    // incremental build to get RNG
    IncrementalBuild IB(pivotLayers2);
    IB.build(false);

    for (int layerIndex = 0; layerIndex < L; layerIndex++) {
        tsl::sparse_set<unsigned int> const& pivots_offline = *((*pivotLayers1)[layerIndex].get_pivotIndices_ptr());
        tsl::sparse_set<unsigned int> const& pivots_incremental = *((*pivotLayers2)[layerIndex].get_pivotIndices_ptr());

        if (pivots_offline != pivots_incremental) {
            delete[] dataPointer;
            BOOST_FAIL("Incremental pivots to not match OfflineBuild");
        }
    }

    std::cout << "  * Done with pivot selection" << std::endl;
    delete[] dataPointer;
    return;
}

void Test_IncrementalBuild::test_familyTree() {
    std::cout << "Unit Test: Family Tree" << std::endl;

    // create dataset
    unsigned int const datasetSize = 800;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0.11395, 0.04995, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension,true);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.build(false);

    std::cout << "L=" <<  L << std::endl;
    // test each coarse layer
    for (int layerIndex = 0; layerIndex < L - 1; layerIndex++) {
        tsl::sparse_set<unsigned int> const& parentIndices = *((*pivotLayers)[layerIndex].get_pivotIndices_ptr());
        float const coarseRadius = effectiveRadiusVector[layerIndex];
        
        // descendants datastructure
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> parentsDescendants_GT{};
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> parentsDescendants_hRNG{};

        // for each layer finer than the current coarse
        for (int finerLayerIndex = layerIndex + 1; finerLayerIndex < L; finerLayerIndex++) {
            tsl::sparse_set<unsigned int> const& childIndices = *((*pivotLayers)[finerLayerIndex].get_pivotIndices_ptr());
            float const fineRadius = effectiveRadiusVector[finerLayerIndex];

            // ancestors structure
            tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> childrensAncestors_GT{};
            tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> childrensAncestors_hRNG{};
            tsl::sparse_set<unsigned int>::const_iterator it0; 
            for (it0 = childIndices.begin(); it0 != childIndices.end(); it0++) {
                unsigned int const childIndex = (*it0);
                childrensAncestors_GT[childIndex] = tsl::sparse_set<unsigned int>{};
                childrensAncestors_hRNG[childIndex] = (*pivotLayers)[finerLayerIndex].get_ancestorPivotIndices(childIndex,layerIndex);
            }

            // iterate through each parent. See if any lower pivots are their children
            tsl::sparse_set<unsigned int>::const_iterator it1; 
            for (it1 = parentIndices.begin(); it1 != parentIndices.end(); it1++) {
                unsigned int const parentIndex = (*it1);
                parentsDescendants_GT[parentIndex] = tsl::sparse_set<unsigned int> {};
                parentsDescendants_hRNG[parentIndex] = (*pivotLayers)[layerIndex].get_descendantPivotIndices(parentIndex,finerLayerIndex);

                // iterate through each lower pivot as potential child
                tsl::sparse_set<unsigned int>::const_iterator it2; 
                for (it2 = childIndices.begin(); it2 != childIndices.end(); it2++) {
                    unsigned int const childIndex = (*it2);
                    float const distance = pivotLayers->getDistance(parentIndex,childIndex);

                    // add parent/child relationship
                    if (distance <= coarseRadius - fineRadius) {
                        parentsDescendants_GT[parentIndex].insert(childIndex);
                        childrensAncestors_GT[childIndex].insert(parentIndex);
                    }
                }
            }

            tsl::sparse_set<unsigned int>::const_iterator it6;
            for (it6 = parentIndices.begin(); it6 != parentIndices.end(); it6++) {
                unsigned int const index = (*it6);
                tsl::sparse_set<unsigned int> const& GT_children = parentsDescendants_GT[index];
                tsl::sparse_set<unsigned int> const& hRNG_children = parentsDescendants_hRNG[index];

                if (GT_children != hRNG_children) {
                    std::cout << "error: layers " << layerIndex << "-" << finerLayerIndex << ". incorrect children on " << index << std::endl;
                    printSet9(GT_children);
                    printSet9(hRNG_children);
                    break;
                }
            }

            if (childrensAncestors_GT != childrensAncestors_hRNG|| parentsDescendants_GT != parentsDescendants_hRNG) {
                std::cout << "Family Tree Error: " << layerIndex << ": " << finerLayerIndex << std::endl;
                delete[] dataPointer;
                BOOST_FAIL("Family Tree does not check out");
            } 
            std::cout << "Passed layer: "<< layerIndex << ": " << finerLayerIndex << std::endl;
        }
    }

    std::cout << "  * Done with family tree" << std::endl;
    delete[] dataPointer;
    return;
}

void Test_IncrementalBuild::test_dmax() {
    std::cout << "Unit Test: Dmax Check" << std::endl;
    bool success = true;

    // create dataset
    unsigned int const datasetSize = 800;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0.11395, 0.04995, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension,true);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.build(false);

    std::cout << "L=" <<  L << std::endl;
    // test each coarse layer
    for (int layerIndex = 0; layerIndex < L - 1; layerIndex++) {
        tsl::sparse_set<unsigned int> const& parentIndices = *((*pivotLayers)[layerIndex].get_pivotIndices_ptr());

        // for each layer finer than the current coarse
        for (int finerLayerIndex = layerIndex + 1; finerLayerIndex < L; finerLayerIndex++) {

            // iterate through each parent. See if any lower pivots are their children
            tsl::sparse_set<unsigned int>::const_iterator it1; 
            for (it1 = parentIndices.begin(); it1 != parentIndices.end(); it1++) {
                unsigned int const parentIndex = (*it1);
                tsl::sparse_set<unsigned int> const& childrenIndices = (*pivotLayers)[layerIndex].get_descendantPivotIndices(parentIndex,finerLayerIndex);
                
                // dmax as max distance to child in finer layer
                float dmax_GT = 0;
                float const dmax_hRNG = (*pivotLayers)[layerIndex].get_dmaxMaxChildDistance(parentIndex,finerLayerIndex);

                // iterate through each child
                tsl::sparse_set<unsigned int>::const_iterator it2;
                for (it2 = childrenIndices.begin(); it2 != childrenIndices.end(); it2++) {
                    unsigned int const childIndex = (*it2);
                    float const distance = pivotLayers->getDistance(parentIndex,childIndex);

                    if (distance > dmax_GT) {
                        dmax_GT = distance;
                    }
                }

                if (dmax_GT != dmax_hRNG) {
                    std::cout << "Dmax Error: " << layerIndex << "-> " << finerLayerIndex << ": " << parentIndex <<  std::endl;
                    std::cout << "GT: " << dmax_GT << ", hRNG: " << dmax_hRNG <<  std::endl;
                    success = false;
                } 
            }
        }
        std::cout << "Passed layer: " << layerIndex  << std::endl;
    }

    if (success == true) {
        std::cout << "dmax test passed!" << std::endl;
    } else {
        std::cout << "dmax test failed!" << std::endl;
        delete[] dataPointer;
        BOOST_FAIL("dmax does not match ground truth!");

    }

    std::cout << "Done with dmax check" << std::endl;
    delete[] dataPointer;
    return;
}

void Test_IncrementalBuild::test_umax() {
    std::cout << "Unit Test: Umax Check" << std::endl;
    bool success = true;

    // create dataset
    unsigned int const datasetSize = 800;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0.11395, 0.04995, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension,true);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.build(false);

    std::cout << "L=" <<  L << std::endl;
    // test each fine layer
    for (int layerIndex = 1; layerIndex < L; layerIndex++) {
        tsl::sparse_set<unsigned int> const& pivotIndices = *((*pivotLayers)[layerIndex].get_pivotIndices_ptr());
        tsl::sparse_map<unsigned int, IndexPairDistanceStruct> const& umax_hRNG = *((*pivotLayers)[layerIndex].get_umaxMaxLinkDistance_ptr());
        tsl::sparse_map<unsigned int, IndexPairDistanceStruct> umax_GT{};
        float const radius = effectiveRadiusVector[layerIndex];

        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
            unsigned int const pivotIndex = (*it1);
            IndexPairDistanceStruct umax_temp = IndexPairDistanceStruct();

            // iterate through each neighbor of pivotINdex to find max link distance
            tsl::sparse_set<unsigned int> const& neighborsList = (*pivotLayers)[layerIndex].get_pivotLayerNeighbors(pivotIndex);
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = neighborsList.begin(); it2 != neighborsList.end(); it2++) {
                unsigned int const neighborIndex = (*it2);
                float const linkDistance = pivotLayers->getDistance(pivotIndex,neighborIndex) - 3*radius;

                if (umax_temp.isInvalid()) {
                    umax_temp = IndexPairDistanceStruct(pivotIndex, neighborIndex, linkDistance);
                } else if (linkDistance > umax_temp.distance) {
                    umax_temp = IndexPairDistanceStruct(pivotIndex, neighborIndex, linkDistance);
                }
            }

            umax_GT[pivotIndex] = umax_temp;
        }

        if (umax_hRNG != umax_GT) {
            std::cout << "umax does not match for layer: " << layerIndex << std::endl;
            delete[] dataPointer;
            BOOST_FAIL("umax does not check out!");
        }

        // now, go through coarser layers and ensure descendant umax is working
        tsl::sparse_map<unsigned int,IndexPairDistanceStruct> umaxD_GT;
        tsl::sparse_map<unsigned int,IndexPairDistanceStruct> umaxD_hRNG;
        for (int coarserLayerIndex = 0; coarserLayerIndex < layerIndex; coarserLayerIndex++) {
            tsl::sparse_set<unsigned int> const& parentIndices = *((*pivotLayers)[coarserLayerIndex].get_pivotIndices_ptr());

            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = parentIndices.begin(); it3 != parentIndices.end(); it3++) {
                unsigned int const parentIndex = (*it3);
                tsl::sparse_set<unsigned int> const& childIndices = (*pivotLayers)[coarserLayerIndex].get_descendantPivotIndices(parentIndex,layerIndex);

                IndexPairDistanceStruct umax_temp = IndexPairDistanceStruct();
                umaxD_hRNG[parentIndex] = (*pivotLayers)[coarserLayerIndex].get_umaxDescendantMaxLinkDistance(parentIndex,layerIndex);

                // check all children to get max link distance 
                tsl::sparse_set<unsigned int>::const_iterator it4;
                for (it4 = childIndices.begin(); it4 != childIndices.end(); it4++) {
                    unsigned int const childIndex = (*it4);
                    IndexPairDistanceStruct const& umax_child = (*pivotLayers)[layerIndex].get_umaxMaxLinkDistance(childIndex);
                    float const linkDistance = pivotLayers->getDistance(parentIndex,childIndex) + umax_child.distance;

                    if (umax_temp.isInvalid()) {
                        umax_temp = IndexPairDistanceStruct(umax_child.index1, umax_child.index2, linkDistance);
                    } else if (linkDistance > umax_temp.distance) {
                        umax_temp = IndexPairDistanceStruct(umax_child.index1, umax_child.index2, linkDistance);
                    }
                }
                umaxD_GT[parentIndex] = umax_temp;
            }

            if (umaxD_hRNG != umaxD_GT) {
                std::cout << "umax descendant does not match for layers " << coarserLayerIndex << ": " << layerIndex << std::endl;

                 tsl::sparse_set<unsigned int>::const_iterator it3;
                for (it3 = parentIndices.begin(); it3 != parentIndices.end(); it3++) {
                    unsigned int const parentIndex = (*it3);
                    if (umaxD_hRNG[parentIndex] != umaxD_GT[parentIndex]) {
                        std::cout << parentIndex << std::endl;
                        std::cout << umaxD_hRNG[parentIndex].index1 << "->" << umaxD_GT[parentIndex].index1 << std::endl;
                        std::cout << umaxD_hRNG[parentIndex].index2 << "->" << umaxD_GT[parentIndex].index2 << std::endl;
                        std::cout << umaxD_hRNG[parentIndex].distance << "->" << umaxD_GT[parentIndex].distance << std::endl;
                    }

                }
                
                delete[] dataPointer;
                BOOST_FAIL("umax descendants does not check out!");
            }
        }

        std::cout << "Passed layer: " << layerIndex  << std::endl;
    }

    std::cout << "Done with Umax Check" << std::endl;
    delete[] dataPointer;
    return;
}

void naiveBerkGraph(float* dataPointer, unsigned int const datasetSize, unsigned int const dimension,
                    tsl::sparse_set<unsigned int> const& indices, float const radius,
                    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& RNG);

void Test_IncrementalBuild::test_neighbors() {
    std::cout << "Unit Test: Neighbors" << std::endl;

    // create dataset
    unsigned int const datasetSize = 1600; //3200;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);

    // 3L hierarchy, 2D, N=void test_buildAndDataStructures();
    // std::vector<float> pivotRadiusVector = std::vector<float>{0.21872, 0};
    // std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0};
    std::vector<float> pivotRadiusVector = std::vector<float>{0.21872, 0.08995, 0.02123, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension,true);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.setNumThreads(4);
    IB.build(true);

    for (int layerIndex = 0; layerIndex < L; layerIndex++) {
        float const radius = effectiveRadiusVector[layerIndex];
        tsl::sparse_set<unsigned int> const& indices = *((*pivotLayers)[layerIndex].get_pivotIndices_ptr());
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& hRNG_RNG = *((*pivotLayers)[layerIndex].get_pivotLayerNeighbors_ptr());

        std::cout << "Comparing hRNG to Ground Truth" << std::endl;

        // get RNG by brute force
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> RNG_GT;
        naiveBerkGraph(dataPointer, datasetSize, dimension, indices, radius, RNG_GT);

        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = indices.begin(); it1 != indices.end(); it1++) {
            unsigned int const index = (*it1);
            tsl::sparse_set<unsigned int> const& GT_neighbors = RNG_GT[index];
            tsl::sparse_set<unsigned int> const& hRNG_neighbors = hRNG_RNG.at(index);

            if (hRNG_neighbors != GT_neighbors) {
                std::cout << "incorrect neighbors on " << index << std::endl;
                std::cout << "GT: "; printSet9(GT_neighbors);
                std::cout << "hRNG: "; printSet9(hRNG_neighbors);
                break;
            }
        }

        // are they equal?
        if (hRNG_RNG != RNG_GT) {
            std::cout << "BG of Layer: " << layerIndex << " does not match GT" << std::endl;
            // printMapSet(hRNG_RNG);
            // printMapSet(RNG_GT);
            delete[] dataPointer;
            BOOST_FAIL("RNG does not match ground truth");
        } else {
            std::cout << "BG of Layer: " << layerIndex << " matches GT!" << std::endl;
        }
    }

    delete[] dataPointer;
    return;
}



void naiveBerkGraph(float* dataPointer, unsigned int const datasetSize, unsigned int const dimension,
                    tsl::sparse_set<unsigned int> const& indices, float const radius,
                    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& RNG) {
    SparseMatrix sparseMatrix(dataPointer, datasetSize, dimension);

    // create vector to sort
    std::vector<unsigned int> indicesVector(indices.begin(), indices.end());
    std::sort(indicesVector.begin(), indicesVector.end());

    RNG.clear();
    std::vector<unsigned int>::const_iterator it1, it2, it3;
    for (it1 = indicesVector.begin(); it1 != indicesVector.end(); it1++) {
        unsigned int index = (*it1);
        RNG[*it1] = tsl::sparse_set<unsigned int>{};
    }

    // consider each index1
    for (it1 = indicesVector.begin(); it1 != indicesVector.end(); it1++) {
        unsigned int const index1 = (*it1);

        // consider each index2 as a neighbor of index1
        for (it2 = indicesVector.begin(); it2 != indicesVector.end(); it2++) {
            unsigned int const index2 = (*it2);
            if (index2 == index1) continue;

            bool flag_isNeighbor = true;
            float const distance12 = sparseMatrix.getDistance(index1, index2);

            // if no lune, no interference --> automatic neighbor
            if (distance12 - 3 * radius <= 0) {
                RNG[index1].insert(index2);
                RNG[index2].insert(index1);
                continue;
            }

            // check for interference by other pivots
            for (it3 = indicesVector.begin(); it3 != indicesVector.end(); it3++) {
                unsigned int const index3 = (*it3);
                if (index3 == index1 || index3 == index2) continue;

                // check both sides for interference
                float const distance13 = sparseMatrix.getDistance(index1, index3);
                if (distance13 < distance12 - 3 * radius) {
                    float const distance23 = sparseMatrix.getDistance(index2, index3);
                    if (distance23 < distance12 - 3 * radius) {
                        flag_isNeighbor = false;
                        break;
                    }
                }
            }

            if (flag_isNeighbor) {
                RNG[index1].insert(index2);
                RNG[index2].insert(index1);
                continue;
            }
        }
    }

    return;
}


void naiveCoarseNeighborsGraph(float* dataPointer, unsigned int const datasetSize, unsigned int const dimension,
                    tsl::sparse_set<unsigned int> const& fineIndices, tsl::sparse_set<unsigned int> const& coarseIndices, float const fineRadius, float const coarseRadius, 
                    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& VGRNG);

void Test_IncrementalBuild::test_coarseNeighbors() {
    std::cout << "Unit Test: CoarseNeighbors" << std::endl;

    // create dataset
    unsigned int const datasetSize = 800;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);

    // 3L hierarchy, 2D, N=void test_buildAndDataStructures();
    //std::vector<float> pivotRadiusVector = std::vector<float>{0.21872, 0};
    //std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0};
    std::vector<float> pivotRadiusVector = std::vector<float>{0.21872, 0.08995, 0.02123, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.build(true);

    for (int layerIndex = 0; layerIndex < L-1; layerIndex++) {
        float const fineRadius = effectiveRadiusVector[layerIndex+1];
        float const coarseRadius = effectiveRadiusVector[layerIndex];
        tsl::sparse_set<unsigned int> const& fineIndices = *((*pivotLayers)[layerIndex+1].get_pivotIndices_ptr());
        tsl::sparse_set<unsigned int> const& coarseIndices = *((*pivotLayers)[layerIndex].get_pivotIndices_ptr());
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> const& hRNG_VGRNG = *((*pivotLayers)[layerIndex+1].get_coarsePivotLayerNeighbors_ptr());

        std::cout << "Comparing hRNG Coarse Neighbors to Ground Truth" << std::endl;

        // get coarse neighbors by brute force
        tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>> VGRNG_GT{};
        naiveCoarseNeighborsGraph(dataPointer, datasetSize, dimension, fineIndices, coarseIndices, fineRadius, coarseRadius, VGRNG_GT);

        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = fineIndices.begin(); it1 != fineIndices.end(); it1++) {
            unsigned int const index = (*it1);
            tsl::sparse_set<unsigned int> const& GT_neighbors = VGRNG_GT[index];
            tsl::sparse_set<unsigned int> const& hRNG_neighbors = hRNG_VGRNG.at(index);

            if (hRNG_neighbors != GT_neighbors) {
                std::cout << "incorrect neighbors on " << index << std::endl;
                printSet9(GT_neighbors);
                printSet9(hRNG_neighbors);
                break;
            }
        }


        // are they equal?
        if (hRNG_VGRNG != VGRNG_GT) {
            std::cout << "V_BG of Layer: " << layerIndex+1 << " does not match GT" << std::endl;
            //printMapSet(hRNG_RNG);
            //printMapSet(RNG_GT);
            delete[] dataPointer;
            BOOST_FAIL("BG does not match ground truth");
        } else {
            std::cout << "BG of Layer: " << layerIndex+1 << " matches GT!" << std::endl;
        }
    }

    delete[] dataPointer;
    return;
}


void naiveCoarseNeighborsGraph(float* dataPointer, unsigned int const datasetSize, unsigned int const dimension,
                    tsl::sparse_set<unsigned int> const& fineIndices, tsl::sparse_set<unsigned int> const& coarseIndices, float const fineRadius, float const coarseRadius, 
                    tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>& VGRNG) {
    SparseMatrix sparseMatrix(dataPointer, datasetSize, dimension);

    VGRNG.clear();
    tsl::sparse_set<unsigned int>::const_iterator it1,it2,it3;
    for (it1 = fineIndices.begin(); it1 != fineIndices.end(); it1++) {
        unsigned int index = (*it1);
        VGRNG[*it1] = tsl::sparse_set<unsigned int>{};
    }

    float radius1 = fineRadius;
    float radius2 = coarseRadius;

    // consider each fine pivot 
    for (it1 = fineIndices.begin(); it1 != fineIndices.end(); it1++) {
        unsigned int const index1 = (*it1);

        // consider each coarse index as a virtual neighbor of index1
        for (it2 = coarseIndices.begin(); it2 != coarseIndices.end(); it2++) {
            unsigned int const index2 = (*it2);
            float const distance12 = sparseMatrix.getDistance(index1, index2);
            bool flag_isNeighbor = true;

            // if no lune, no interference --> automatic virtual neighbor
            if (distance12 - 2*radius1 - radius2 <= 0 || distance12 - radius1 - 2*radius2 <= 0) {
                VGRNG[index1].insert(index2);
                continue;
            }

            // check for interference by other pivots
            for (it3 = coarseIndices.begin(); it3 != coarseIndices.end(); it3++) {
                unsigned int const index3 = (*it3);
                if (index3 == index2) continue;

                // check both sides for interference
                float const distance13 = sparseMatrix.getDistance(index1, index3);
                if (distance13 < distance12 - 2*radius1 - radius2) {
                    float const distance23 = sparseMatrix.getDistance(index2, index3);
                    if (distance23 < distance12 - radius1 - 2*radius2) {
                        flag_isNeighbor = false;
                        break;
                    }
                }
            }

            if (flag_isNeighbor) {
                VGRNG[index1].insert(index2);
            }
        }
    }

    return;
}




void Test_IncrementalBuild::test_vmax() {
    std::cout << "Unit Test: Vmax Check" << std::endl;
    bool success = true;

    // create dataset
    unsigned int const datasetSize = 800;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::vector<float> pivotRadiusVector = std::vector<float>{0.46079, 0.21872, 0.11395, 0.04995, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);
    int L = effectiveRadiusVector.size();

    // create pivot layer hierarchy, initializes sparse matrix
    std::shared_ptr<PivotLayerHierarchy> pivotLayers =
        std::make_shared<PivotLayerHierarchy>(effectiveRadiusVector, "Berk", dataPointer, datasetSize, dimension, true);

    // offline build to get RNG
    IncrementalBuild IB(pivotLayers);
    IB.build(false);

    std::cout << "L=" <<  L << std::endl;
    // test each fine layer
    for (int layerIndex = 1; layerIndex < L; layerIndex++) {
        int const coarserLayerIndex = layerIndex - 1;
        tsl::sparse_set<unsigned int> const& pivotIndices = *((*pivotLayers)[layerIndex].get_pivotIndices_ptr());
        tsl::sparse_map<unsigned int, IndexPairDistanceStruct> const& vmax_hRNG = *((*pivotLayers)[layerIndex].get_vmaxMaxCoarseLinkDistance_ptr());
        tsl::sparse_map<unsigned int, IndexPairDistanceStruct> vmax_GT{};
        float const fineRadius = effectiveRadiusVector[layerIndex];
        float const coarseRadius = effectiveRadiusVector[coarserLayerIndex];
        
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
            unsigned int const pivotIndex = (*it1);
            IndexPairDistanceStruct vmax_temp;

            // iterate through each neighbor of pivotINdex to find max link distance
            tsl::sparse_set<unsigned int> const& coarseNeighborsList = (*pivotLayers)[layerIndex].get_coarsePivotLayerNeighbors(pivotIndex);
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = coarseNeighborsList.begin(); it2 != coarseNeighborsList.end(); it2++) {
                unsigned int const coarseIndex = (*it2);
                float const distance = pivotLayers->getDistance(pivotIndex,coarseIndex);
                float const linkDistance = pivotLayers->_luneRadius(distance,fineRadius,coarseRadius);

                if (vmax_temp.isInvalid()) {
                    vmax_temp = IndexPairDistanceStruct(pivotIndex, coarseIndex, linkDistance);
                } else if (linkDistance > vmax_temp.distance) {
                    vmax_temp = IndexPairDistanceStruct(pivotIndex, coarseIndex, linkDistance);
                }
            }

            vmax_GT[pivotIndex] = vmax_temp;
        }

        if (vmax_hRNG != vmax_GT) {
            std::cout << "vmax does not match for layer: " << layerIndex << std::endl;
            delete[] dataPointer;
            BOOST_FAIL("vmax does not check out!");
        }

        // now, go through coarser layers and ensure descendant vmax is working
        tsl::sparse_map<unsigned int,IndexPairDistanceStruct> vmaxD_GT;
        tsl::sparse_map<unsigned int,IndexPairDistanceStruct> vmaxD_hRNG;
        for (int coarserLayerIndex = 0; coarserLayerIndex < layerIndex; coarserLayerIndex++) {
            tsl::sparse_set<unsigned int> const& parentIndices = *((*pivotLayers)[coarserLayerIndex].get_pivotIndices_ptr());

            tsl::sparse_set<unsigned int>::const_iterator it3;
            for (it3 = parentIndices.begin(); it3 != parentIndices.end(); it3++) {
                unsigned int const parentIndex = (*it3);
                tsl::sparse_set<unsigned int> const& childIndices = (*pivotLayers)[coarserLayerIndex].get_descendantPivotIndices(parentIndex,layerIndex);

                IndexPairDistanceStruct vmax_temp = IndexPairDistanceStruct();
                vmaxD_hRNG[parentIndex] = (*pivotLayers)[coarserLayerIndex].get_vmaxDescendantMaxCoarseLinkDistance(parentIndex,layerIndex);

                // check all children to get max link distance 
                tsl::sparse_set<unsigned int>::const_iterator it4;
                for (it4 = childIndices.begin(); it4 != childIndices.end(); it4++) {
                    unsigned int const childIndex = (*it4);
                    IndexPairDistanceStruct const& vmax_child = (*pivotLayers)[layerIndex].get_vmaxMaxCoarseLinkDistance(childIndex);
                    float const linkDistance = pivotLayers->getDistance(parentIndex,childIndex) + vmax_child.distance;

                    if (vmax_temp.isInvalid()) {
                        vmax_temp = IndexPairDistanceStruct(vmax_child.index1, vmax_child.index2, linkDistance);
                    } else if (linkDistance > vmax_temp.distance) {
                        vmax_temp = IndexPairDistanceStruct(vmax_child.index1, vmax_child.index2, linkDistance);
                    }
                }
                vmaxD_GT[parentIndex] = vmax_temp;
            }

            if (vmaxD_hRNG != vmaxD_GT) {
                std::cout << "vmax descendant does not match for layers " << coarserLayerIndex << ": " << layerIndex << std::endl;
                delete[] dataPointer;
                BOOST_FAIL("vmax descendants does not check out!");
            }
        }

        std::cout << "Passed layer: " << layerIndex  << std::endl;
    }

    std::cout << "Done with Vmax Check" << std::endl;
    delete[] dataPointer;
    return;
}




