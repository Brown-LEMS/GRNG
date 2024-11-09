#define BOOST_TEST_DYN_LINK  // deals with no main() issue
#define BOOST_TEST_MODULE PivotSelectionTest
#include <boost/test/unit_test.hpp>

#include "constant-radius-pivot-selection.hpp"
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

BOOST_AUTO_TEST_SUITE(pivot_selection)
BOOST_AUTO_TEST_CASE(correctness) {
    unsigned int const datasetSize = 6400;
    unsigned int const dimension = 2;
    float* dataPointer;
    Datasets::uniformDataset(dataPointer, dimension, datasetSize, 3);
    std::shared_ptr<SparseMatrix> SM = std::make_shared<SparseMatrix>(dataPointer, datasetSize, dimension);

    // 5L hierarchy, 2D, N=12800
    std::vector<float> pivotRadiusVector = std::vector<float>{0.22006, 0.12608, 0.07401, 0.03783, 0};
    std::vector<float> effectiveRadiusVector = getEffectiveRadiusVector(pivotRadiusVector);

    // constant radius pivot selection
    ConstantRadiusPivotSelection ps(datasetSize, effectiveRadiusVector, SM, false);
    std::vector<PivotSelectionStruct> const selectionStruct = ps.selectPivots();

    BOOST_CHECK(selectionStruct[0].pivotIndices->size() == 47);
    BOOST_CHECK(selectionStruct[1].pivotIndices->size() == 131);
    BOOST_CHECK(selectionStruct[2].pivotIndices->size() == 387);
    BOOST_CHECK(selectionStruct[3].pivotIndices->size() == 1375);
    BOOST_CHECK(selectionStruct[4].pivotIndices->size() == datasetSize);

    // check that each pivot is a child of a pivot in the above layer...
    for (int layerIndex = 1; layerIndex < pivotRadiusVector.size(); layerIndex++) {
        tsl::sparse_set<unsigned int> const& indices = *selectionStruct[layerIndex].pivotIndices;
        float const pivotRadius = effectiveRadiusVector[layerIndex];

        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = indices.begin(); it1 != indices.end(); it1++) {
            unsigned int const pivotIndex = (*it1);
            bool orphan = true;

            tsl::sparse_set<unsigned int> const& coarseIndices = *selectionStruct[layerIndex - 1].pivotIndices;
            float const coarseRadius = effectiveRadiusVector[layerIndex - 1];

            // check to make sure each pivot has a parent in its above layer
            tsl::sparse_set<unsigned int>::const_iterator it2;
            for (it2 = coarseIndices.begin(); it2 != coarseIndices.end(); it2++) {
                unsigned int const coarseIndex = (*it2);
                float const distance = SM->getDistance(pivotIndex, coarseIndex);
                if (coarseRadius > distance + pivotRadius) {
                    orphan = false;
                    break;
                }
            }
            BOOST_REQUIRE(!orphan);
        }
    }

    delete[] dataPointer;
    return;
}

BOOST_AUTO_TEST_SUITE_END()
