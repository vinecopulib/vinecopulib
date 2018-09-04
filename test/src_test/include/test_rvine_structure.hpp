// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/vinecop/rvine_structure.hpp>


namespace test_rvine_structure {

using namespace vinecopulib;

TEST(rvine_structure, triangular_array_works) {

    TriangularArray<size_t> my_rvm(7);
    my_rvm[0] = {2, 1, 3, 4, 6, 5};
    my_rvm[1] = {3, 1, 2, 4, 5};
    my_rvm[2] = {1, 4, 3, 2};
    my_rvm[3] = {1, 3, 2};
    my_rvm[4] = {1, 2};
    my_rvm[5] = {1};

    EXPECT_NO_THROW(my_rvm.str());

    std::ostringstream oss;
    oss << my_rvm;
    EXPECT_EQ(oss.str(), my_rvm.str());
    
    std::vector<size_t> myvec = {1, 2};
    EXPECT_NO_THROW(my_rvm.set_column(4, myvec));
    EXPECT_ANY_THROW(my_rvm.set_column(6, myvec));
    EXPECT_ANY_THROW(my_rvm.set_column(3, myvec));
        
    my_rvm.truncate(5);
    EXPECT_EQ(my_rvm.get_trunc_lvl(), 5);
    EXPECT_EQ(my_rvm[0].size(), 5);
    
    my_rvm.reduce(4);
    EXPECT_EQ(my_rvm.get_trunc_lvl(), 3);
    EXPECT_EQ(my_rvm[0].size(), 3);
    
    my_rvm.truncate(1);
    my_rvm.reduce(3);
    EXPECT_EQ(my_rvm.get_trunc_lvl(), 1);
    EXPECT_EQ(my_rvm[0].size(), 1);
    EXPECT_EQ(my_rvm[1].size(), 1);
}

TEST(rvine_structure, can_convert_to_natural_order) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;

    TriangularArray<size_t> true_no_array(7);
    true_no_array[0] = {2, 1, 3, 4, 6, 5};
    true_no_array[1] = {3, 1, 2, 4, 5};
    true_no_array[2] = {1, 4, 3, 2};
    true_no_array[3] = {1, 3, 2};
    true_no_array[4] = {1, 2};
    true_no_array[5] = {1};

    RVineStructure rvine_structure(mat);
    EXPECT_EQ(rvine_structure.get_struct_array(), true_no_array);
}

TEST(rvine_structure, max_array_is_correct) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;

    TriangularArray<size_t> true_max_array(7);
    true_max_array[0] = {2, 2, 3, 4, 6, 6};
    true_max_array[1] = {3, 3, 3, 4, 5};
    true_max_array[2] = {1, 4, 4, 4};
    true_max_array[3] = {1, 3, 3};
    true_max_array[4] = {1, 2};
    true_max_array[5] = {1};

    RVineStructure rvine_structure(mat);
    EXPECT_EQ(rvine_structure.get_max_array(), true_max_array);
    
    rvine_structure.truncate(2);
    EXPECT_EQ(rvine_structure.get_trunc_lvl(), 2);
    EXPECT_EQ(rvine_structure.get_max_array()[0].size(), 2);
}

TEST(rvine_structure, needed_hfunc1_is_correct) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;

    TriangularArray<size_t> true_hfunc1(7);
    true_hfunc1[0] = {0, 0, 0, 0, 0, 0};
    true_hfunc1[1] = {0, 0, 0, 0, 1};
    true_hfunc1[2] = {0, 0, 0, 0};
    true_hfunc1[3] = {0, 1, 1};
    true_hfunc1[4] = {1, 1};
    true_hfunc1[5] = {1};

    RVineStructure rvine_structure(mat);
    EXPECT_EQ(rvine_structure.get_needed_hfunc1(), true_hfunc1);
    
    rvine_structure.truncate(2);
    EXPECT_EQ(rvine_structure.get_needed_hfunc1()[0].size(), 2);
}

TEST(rvine_structure, needed_hfunc2_is_correct) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;

    TriangularArray<size_t> true_hfunc2(7);
    true_hfunc2[0] = {1, 1, 1, 1, 1, 0};
    true_hfunc2[1] = {1, 1, 1, 1, 0};
    true_hfunc2[2] = {1, 1, 1, 1};
    true_hfunc2[3] = {1, 1, 1};
    true_hfunc2[4] = {1, 1};
    true_hfunc2[5] = {1};

    RVineStructure rvine_structure(mat);
    EXPECT_EQ(rvine_structure.get_needed_hfunc2(), true_hfunc2);
    
    rvine_structure.truncate(2);
    EXPECT_EQ(rvine_structure.get_needed_hfunc2()[0].size(), 2);
}

TEST(rvine_structure, construct_d_vine_struct_is_correct) {

    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
        true_d_vine_matrix(7, 7);
    true_d_vine_matrix << 4, 1, 5, 3, 2, 7, 7,
        1, 5, 3, 2, 7, 2, 0,
        5, 3, 2, 7, 3, 0, 0,
        3, 2, 7, 5, 0, 0, 0,
        2, 7, 1, 0, 0, 0, 0,
        7, 4, 0, 0, 0, 0, 0,
        6, 0, 0, 0, 0, 0, 0;

    std::vector<size_t> order = {7, 2, 3, 5, 1, 4, 6};
    RVineStructure rvine_structure(order);
    EXPECT_EQ(rvine_structure.get_matrix(), true_d_vine_matrix);
}

TEST(rvine_structure, rvine_struct_sanity_checks_work) {
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
    mat << 5, 2, 6, 6, 6, 6, 6,
        6, 6, 1, 2, 5, 5, 0,
        2, 5, 2, 5, 2, 0, 0,
        1, 1, 5, 1, 0, 0, 0,
        3, 7, 7, 0, 0, 0, 0,
        7, 3, 0, 0, 0, 0, 0,
        4, 0, 0, 0, 0, 0, 0;

    // should pass without errors
    auto rvm = RVineStructure(mat);
    auto wrong_mat = mat;

    // must be quadratic
    wrong_mat = mat.block(0, 0, 4, 5);
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

    // lower right triangle must contain zeros
    wrong_mat = mat;
    wrong_mat(6, 6) = 1;
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

    // upper left triangle must only contain 1, ..., d
    wrong_mat = mat;
    wrong_mat(0, 0) = 9;
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));
    wrong_mat(0, 0) = 0;
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

    // diagonal elements cannot appear further to the right
    wrong_mat = mat;
    wrong_mat(0, 1) = 4;
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

    // all numbers in a column most appear in each column to the left
    wrong_mat = mat;
    wrong_mat(0, 0) = 4;
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

    // proximity condition
    wrong_mat = mat;
    wrong_mat(3, 1) = 7;
    wrong_mat(4, 1) = 1;
    EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));
}
}
