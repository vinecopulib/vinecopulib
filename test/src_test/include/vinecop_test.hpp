/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <RInside.h>
#include <RcppEigen.h>

class VinecopTest {
public:
    VinecopTest() {
        R.parseEval("library(VineCopula)");
        R.parseEval("set.seed(5)");
        R.parseEval("mat <- matrix(c(4, 0, 0, 0, 0, 0, 0,"
                                   " 7, 3, 0, 0, 0, 0, 0,"
                                   " 3, 7, 7, 0, 0, 0, 0,"
                                   " 1, 1, 5, 1, 0, 0, 0,"
                                   " 2, 5, 2, 5, 2, 0, 0,"
                                   " 6, 6, 1, 2, 5, 5, 0,"
                                   " 5, 2, 6, 6, 6, 6, 6),"
                                   " 7, 7, byrow = TRUE)");
        R.parseEval("fam <- par <- matrix(0, 7, 7)");
        R.parseEval("fam[lower.tri(fam)] <- 23");
        R.parseEval("par[lower.tri(par)] <- -3");
        R.parseEval("model <- RVineMatrix(mat, fam, par)");
        R.parseEval("u <- RVineSim(1000, model)");
        R.parseEval("fit <- RVineStructureSelect(u, familyset = 3)");
        u = Rcpp::as<MatXd>(R.parseEval("u"));
        // C++ representation reverses rows
        matrix = Rcpp::as<MatXi>(R.parseEval("mat")).colwise().reverse();
    }
    
    RInside R;
    MatXd u;
    MatXi matrix;
};


