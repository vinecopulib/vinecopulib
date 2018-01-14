// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/vinecop/rvine_matrix.hpp>
#include <iostream>

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
    EXPECT_EQ(rvine_matrix.get_natural_order(), true_no_matrix);
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
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> true_d_vine_matrix(7,
                                                                             7);
    true_d_vine_matrix << 4, 1, 5, 3, 2, 7, 7,
        1, 5, 3, 2, 7, 2, 0,
        5, 3, 2, 7, 3, 0, 0,
        3, 2, 7, 5, 0, 0, 0,
        2, 7, 1, 0, 0, 0, 0,
        7, 4, 0, 0, 0, 0, 0,
        6, 0, 0, 0, 0, 0, 0;
    EXPECT_EQ(RVineMatrix::construct_d_vine_matrix(order), true_d_vine_matrix);
}

TEST(rvine_matrix, belongs_to_structure_is_correct) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;
    RVineMatrix rvine_matrix(mat);

    // conditioned set size is equal to 2
    std::vector<size_t> conditioned0 = {0, 1, 2, 3};
    std::vector<size_t> conditioning0 = {0, 1, 2, 3};
    EXPECT_ANY_THROW(rvine_matrix.belongs_to_structure(conditioned0,
                                                       conditioning0));

    std::vector<size_t> conditioned1 = {4, 5};
    std::vector<size_t> conditioning1(0);
    ASSERT_TRUE(rvine_matrix.belongs_to_structure(conditioned1,
                                                  conditioning1));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned1,
                                                   conditioning0));
    conditioned1[0] = 5;
    conditioned1[1] = 4;
    // only allow for correctly ordered conditioned set
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned1,
                                                   conditioning1));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned1,
                                                   conditioning0));
    conditioned1[0] = 5;
    conditioned1[1] = 6;
    ASSERT_TRUE(rvine_matrix.belongs_to_structure(conditioned1,
                                                  conditioning1));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned1,
                                                   conditioning0));

    std::vector<size_t> conditioned2 = {4, 6};
    std::vector<size_t> conditioning2 = {5};
    ASSERT_TRUE(rvine_matrix.belongs_to_structure(conditioned2,
                                                  conditioning2));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned2,
                                                   conditioning1));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned2,
                                                   conditioning0));
    conditioning2[0] = 6;
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned2,
                                                   conditioning2));
    conditioned2[0] = 2;
    conditioned2[1] = 5;
    ASSERT_TRUE(rvine_matrix.belongs_to_structure(conditioned2,
                                                  conditioning2));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned2,
                                                   conditioning1));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned2,
                                                   conditioning0));

    std::vector<size_t> conditioned3 = {1, 5};
    std::vector<size_t> conditioning3 = {6, 2};
    ASSERT_TRUE(rvine_matrix.belongs_to_structure(conditioned3,
                                                  conditioning3));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned3,
                                                   conditioning1));
    ASSERT_FALSE(rvine_matrix.belongs_to_structure(conditioned3,
                                                   conditioning0));

}

TEST(rvine_matrix, rvine_matrix_sanity_checks_work) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;
    // should pass without errors
    auto rvm = RVineMatrix(mat);

    // must be quadratic
    auto wrong_mat = mat;
    wrong_mat = mat.block(0, 0, 4, 5);
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));

    // lower right triangle must contain zeros
    wrong_mat = mat;
    wrong_mat(6, 6) = 1;
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));

    // upper left triangle must only contain 1, ..., d
    wrong_mat = mat;
    wrong_mat(0, 0) = 9;
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));
    wrong_mat(0, 0) = 0;
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));

    // diagonal elements cannot appear further to the right
    wrong_mat = mat;
    wrong_mat(0, 1) = 4;
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));

    // all numbers in a column most appear in each column to the left
    wrong_mat = mat;
    wrong_mat(0, 0) = 4;
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));

    // proximity condition
    wrong_mat = mat;
    wrong_mat(3, 1) = 7;
    wrong_mat(4, 1) = 1;
    EXPECT_ANY_THROW(rvm = RVineMatrix(wrong_mat));

    // row and col should be smaller than d_
    EXPECT_ANY_THROW(rvm.get_element(8, 0));
    EXPECT_ANY_THROW(rvm.get_element(0, 8));
}
}
