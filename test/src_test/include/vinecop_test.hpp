/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VINECOPLIB_VINECOP_TEST_HPP
#define VINECOPLIB_VINECOP_TEST_HPP

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


#endif
