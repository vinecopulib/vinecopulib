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
#include "include/rvine_matrix.hpp"

namespace {
    TEST(rvine_matrix, can_convert_to_natural_order) {
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        MatXi true_no_matrix(7, 7);
        true_no_matrix << 7, 0, 0, 0, 0, 0, 0,
                          5, 6, 0, 0, 0, 0, 0,
                          6, 5, 5, 0, 0, 0, 0,
                          4, 4, 2, 4, 0, 0, 0,
                          3, 2, 3, 2, 3, 0, 0,
                          1, 1, 4, 3, 2, 2, 0,
                          2, 3, 1, 1, 1, 1, 1;
                          
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.in_natural_order(), true_no_matrix);
    }
    
    TEST(rvine_matrix, max_mat_is_correct) {
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        MatXi true_max_matrix(7, 7);
        true_max_matrix << 7, 0, 0, 0, 0, 0, 0,
                           6, 6, 0, 0, 0, 0, 0,
                           6, 5, 5, 0, 0, 0, 0,
                           4, 4, 4, 4, 0, 0, 0,
                           3, 3, 4, 3, 3, 0, 0,
                           2, 3, 4, 3, 2, 2, 0,
                           2, 3, 1, 1, 1, 1, 1;
                           
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.get_max_matrix(), true_max_matrix);
    }
    
    TEST(rvine_matrix, needed_hfunc1_is_correct) {
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        MatXb true_hfunc1(7, 7);
        true_hfunc1 << 0, 0, 0, 0, 0, 0, 0,
                       1, 0, 0, 0, 0, 0, 0,
                       1, 1, 1, 0, 0, 0, 0,
                       1, 1, 1, 1, 0, 0, 0,
                       1, 1, 1, 1, 1, 0, 0,
                       1, 1, 1, 1, 1, 1, 0,
                       1, 1, 1, 1, 1, 1, 0;
                       
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.get_needed_hfunc1(), true_hfunc1);
    }
    
    TEST(rvine_matrix, needed_hfunc2_is_correct) {
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        MatXb true_hfunc2(7, 7);
        true_hfunc2 << 0, 0, 0, 0, 0, 0, 0,
                       0, 1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 1, 0, 0, 0,
                       0, 0, 0, 1, 1, 0, 0,
                       0, 0, 0, 0, 1, 1, 0,
                       0, 0, 0, 0, 0, 0, 0;
                       
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.get_needed_hfunc2(), true_hfunc2);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
