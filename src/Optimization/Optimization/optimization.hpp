#ifndef Optimization_hpp
#define Optimization_hpp
/*
Copyright 2022, Brown University, Providence, RI.

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
// 2022-06-04
#include <string>
#include <vector>

#include "hrng-data.hpp"

class Optimization {
  public:
    Optimization(unsigned int const dimension, float *const &dataPointer, unsigned int const datasetSize, float *const &testPointer,
                 unsigned int const testsetSize, int const numThreads, std::string luneType);
    Optimization(unsigned int const dimension, float *const &dataPointer, unsigned int const datasetSize, float *const &testPointer,
                 unsigned int const testsetSize, int const numThreads, std::string luneType, bool cacheAll);
    ~Optimization(){};

    double static objectiveFunction(const std::vector<double> &x, std::vector<double> &grad, void *hRNG_Data_Void);

    bool localOptimization(std::string method, unsigned n, std::vector<double> &x, double &f, std::vector<double> lowerBounds,
                           std::vector<double> upperBounds, std::vector<double> stepSize, std::vector<double> xTolVector, double fTol, int maxEval,
                           std::string resultsSaverPath);

    hRNGData _data;
};

#endif  // Optimization_hpp