#ifndef PivotSelectionStruct_hpp
#define PivotSelectionStruct_hpp
/*
Copyright 2019, Brown University, Providence, RI.
                        All Rights Reserved
Permission to use, copy, modify, and distribute this software and
its documentation for any purpose other than its incorporation into a
commercial product or service is hereby granted without fee, provided
that the above copyright notice appear in all copies and that both
that copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.
BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
// Cole Foster
// February 1st, 2022

// import libraries
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"

/**
 * @brief Struct to hold pivot indices and radii for a particular layer. Used as return object for pivot selection.
 *
 */
struct PivotSelectionStruct {
  public:
    PivotSelectionStruct(std::shared_ptr<tsl::sparse_set<unsigned int>> const& pivotIndices,
                         std::shared_ptr<tsl::sparse_map<unsigned int, float>> const& pivotRadii) : pivotIndices(pivotIndices), pivotRadii(pivotRadii){};
    ~PivotSelectionStruct(){};

    std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices;
    std::shared_ptr<tsl::sparse_map<unsigned int, float>> pivotRadii;
};

#endif  // PivotSelectionStruct_hpp