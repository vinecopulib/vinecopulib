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
            pair_copulas[i] = Bicop::create(3, VecXd::Constant(1, 3.0), 90);
        }
        MatXd u(10, 7);
        u << 0.72, 0.82, 0.18, 0.39, 0.64, 0.27, 0.56,
             0.57, 0.25, 0.98, 0.54, 0.44, 0.66, 0.28,
             0.29, 0.33, 0.39, 0.51, 0.26, 0.77, 0.16,
             0.89, 0.58, 0.74, 0.53, 0.83, 0.21, 0.63,
             0.13, 0.18, 0.57, 0.85, 0.12, 0.87, 0.14,
             0.76, 0.21, 0.92, 0.30, 0.55, 0.67, 0.18,
             0.26, 0.25, 0.59, 0.56, 0.26, 0.78, 0.19,
             0.20, 0.53, 0.25, 0.52, 0.34, 0.69, 0.54,
             0.57, 0.75, 0.23, 0.38, 0.69, 0.32, 0.90,
             0.82, 0.54, 0.59, 0.55, 0.76, 0.33, 0.56;
        VecXd true_pdf(10);
        true_pdf << 8.214122e-07, 9.749805e+03, 7.600675e-20, 1.225511e+03,
         2.496552e+04, 5.193771e+05, 1.387955e-20, 2.960827e+05, 7.874391e-04,
         6.707183e+00;
        Vinecop vinecop(pair_copulas, mat);
        ASSERT_TRUE(vinecop.pdf(u).isApprox(true_pdf, 1e-4));
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
