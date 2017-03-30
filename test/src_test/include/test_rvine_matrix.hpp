// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/vinecop/rvine_matrix.hpp>

namespace test_rvine_matrix {
    using namespace vinecopulib;

    TEST(rvine_matrix, can_convert_to_natural_order) {
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
                6, 6, 1, 2, 5, 5, 0,
                2, 5, 2, 5, 2, 0, 0,
                1, 1, 5, 1, 0, 0, 0,
                3, 7, 7, 0, 0, 0, 0,
                7, 3, 0, 0, 0, 0, 0,
                4, 0, 0, 0, 0, 0, 0;
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> true_no_matrix(7, 7);
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
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
                6, 6, 1, 2, 5, 5, 0,
                2, 5, 2, 5, 2, 0, 0,
                1, 1, 5, 1, 0, 0, 0,
                3, 7, 7, 0, 0, 0, 0,
                7, 3, 0, 0, 0, 0, 0,
                4, 0, 0, 0, 0, 0, 0;
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> true_max_matrix(7, 7);
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
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
                6, 6, 1, 2, 5, 5, 0,
                2, 5, 2, 5, 2, 0, 0,
                1, 1, 5, 1, 0, 0, 0,
                3, 7, 7, 0, 0, 0, 0,
                7, 3, 0, 0, 0, 0, 0,
                4, 0, 0, 0, 0, 0, 0;
        MatrixXb true_hfunc1(7, 7);
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
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
        mat << 5, 2, 6, 6, 6, 6, 6,
                6, 6, 1, 2, 5, 5, 0,
                2, 5, 2, 5, 2, 0, 0,
                1, 1, 5, 1, 0, 0, 0,
                3, 7, 7, 0, 0, 0, 0,
                7, 3, 0, 0, 0, 0, 0,
                4, 0, 0, 0, 0, 0, 0;
        MatrixXb true_hfunc2(7, 7);
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
        Eigen::Matrix<size_t, Eigen::Dynamic, 1> order(7);
        order << 7, 2, 3, 5, 1, 4, 6;
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> true_d_vine_matrix(7, 7);
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
