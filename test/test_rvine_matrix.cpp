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

#include "gtest/gtest.h"
#include "include/rvine_matrix.hpp"

namespace {
    TEST(rvine_matrix, can_convert_to_natural_order) {
        MatXi mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
               6, 6, 1, 2, 5, 5, 0,
               2, 5, 2, 5, 2, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               4, 0, 0, 0, 0, 0, 0;
        MatXi true_no_matrix(7, 7);
        true_no_matrix << 2, 3, 1, 1, 1, 1, 1,
                          1, 1, 4, 3, 2, 2, 0,
                          3, 2, 3, 2, 3, 0, 0,
                          4, 4, 2, 4, 0, 0, 0,
                          6, 5, 5, 0, 0, 0, 0,
                          5, 6, 0, 0, 0, 0, 0,
                          7, 0, 0, 0, 0, 0, 0;
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.in_natural_order(), true_no_matrix);
    }
    
    TEST(rvine_matrix, max_mat_is_correct) {
        MatXi mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
               6, 6, 1, 2, 5, 5, 0,
               2, 5, 2, 5, 2, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               4, 0, 0, 0, 0, 0, 0;
        MatXi true_max_matrix(7, 7);
        true_max_matrix << 2, 3, 1, 1, 1, 1, 1,
                           2, 3, 4, 3, 2, 2, 0,
                           3, 3, 4, 3, 3, 0, 0,
                           4, 4, 4, 4, 0, 0, 0,
                           6, 5, 5, 0, 0, 0, 0,
                           6, 6, 0, 0, 0, 0, 0,
                           7, 0, 0, 0, 0, 0, 0;
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.get_max_matrix(), true_max_matrix);
    }
    
    TEST(rvine_matrix, needed_hfunc1_is_correct) {
        MatXi mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
               6, 6, 1, 2, 5, 5, 0,
               2, 5, 2, 5, 2, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               4, 0, 0, 0, 0, 0, 0;
        MatXb true_hfunc1(7, 7);
        true_hfunc1 << 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 1, 1, 0,
                       0, 0, 0, 1, 1, 0, 0,
                       0, 0, 0, 1, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0,
                       0, 1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0;
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.get_needed_hfunc1(), true_hfunc1);
    }
    
    TEST(rvine_matrix, needed_hfunc2_is_correct) {
        MatXi mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
               6, 6, 1, 2, 5, 5, 0,
               2, 5, 2, 5, 2, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               4, 0, 0, 0, 0, 0, 0;
        MatXb true_hfunc2(7, 7);
        true_hfunc2 << 1, 1, 1, 1, 1, 1, 0,
                       1, 1, 1, 1, 1, 1, 0,
                       1, 1, 1, 1, 1, 0, 0,
                       1, 1, 1, 1, 0, 0, 0,
                       1, 1, 1, 0, 0, 0, 0,
                       1, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0;   
        RVineMatrix rvine_matrix(mat);
        EXPECT_EQ(rvine_matrix.get_needed_hfunc2(), true_hfunc2);
    }
    
    TEST(rvine_matrix, construct_d_vine_matrix_is_correct) {
        VecXi order(7);
        order << 7, 2, 3, 5, 1, 4, 6;
        MatXi true_d_vine_matrix(7, 7);
        true_d_vine_matrix << 4, 1, 5, 3, 2, 7, 7,
                              1, 5, 3, 2, 7, 2, 0,
                              5, 3, 2, 7, 3, 0, 0,
                              3, 2, 7, 5, 0, 0, 0,
                              2, 7, 1, 0, 0, 0, 0,
                              7, 4, 0, 0, 0, 0, 0,
                              6, 0, 0, 0, 0, 0, 0;
        EXPECT_EQ(RVineMatrix::construct_d_vine_matrix(order), true_d_vine_matrix);
     }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
