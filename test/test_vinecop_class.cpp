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

#include "gtest/gtest.h"
#include "include/vinecop_class.hpp"

namespace {
    TEST(vinecop_class, constructors_without_error) {
        Vinecop vinecop_default(5);
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        std::vector<BicopPtr> pair_copulas(28);
        for (int i = 0; i < 28; ++i) {
            pair_copulas[i] = Bicop::create(1, 0.5);
        }
        Vinecop vinecop_parametrized(pair_copulas, mat);
    }
    
    TEST(vinecop_class, pdf_is_correct) {
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        std::vector<BicopPtr> pair_copulas(28);
        for (int i = 0; i < 28; ++i) {
            pair_copulas[i] = Bicop::create(3, VecXd::Constant(1, 3.0), 0);
        }
        MatXd u(10, 7);
        u << 0.31, 0.28, 0.27, 0.33, 0.26, 0.27, 0.25,
             0.62, 0.56, 0.70, 0.52, 0.74, 0.66, 0.58,
             0.94, 0.83, 0.77, 0.94, 0.74, 0.77, 0.72,
             0.24, 0.21, 0.23, 0.18, 0.29, 0.21, 0.23,
             0.72, 0.72, 0.81, 0.89, 0.70, 0.87, 0.70,
             0.86, 0.66, 0.90, 0.60, 0.85, 0.67, 0.68,
             0.84, 0.80, 0.75, 0.79, 0.77, 0.78, 0.73,
             0.68, 0.85, 0.60, 0.77, 0.68, 0.69, 0.76,
             0.36, 0.43, 0.31, 0.40, 0.36, 0.32, 0.52,
             0.44, 0.38, 0.38, 0.34, 0.45, 0.33, 0.41;
        VecXd true_pdf(10);
        true_pdf << 6488041, 158077, 786, 153678, 119496, 
                    62551, 322698, 1040453, 37891, 23475;
        Vinecop vinecop(pair_copulas, mat);
        ASSERT_TRUE(vinecop.pdf(u).isApprox(true_pdf, 1.0));
        // TODO: This does not work for 90 / 270 rotations because nans appear. 
        // It comes from rounding errors in the hfunctions that have to be fixed.
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
