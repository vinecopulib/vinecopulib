// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "include/bicop_class.hpp"

namespace test_bicop_class {
    using namespace vinecopulib;

    TEST(bicop_class, creates_right_copula) {
        BicopPtr bicop = Bicop::create(1, Eigen::VectorXd::Ones(1), 90);
        EXPECT_EQ(bicop->get_family(), 1);
        EXPECT_EQ(bicop->get_rotation(), 90);
        EXPECT_EQ(bicop->get_parameters(), Eigen::VectorXd::Ones(1));
    }

    TEST(bicop_class, catches_wrong_parameter_size) {
        EXPECT_ANY_THROW(Bicop::create(0, Eigen::VectorXd::Zero(1), 0));
        EXPECT_ANY_THROW(Bicop::create(1, Eigen::VectorXd::Zero(0), 0));
        EXPECT_ANY_THROW(Bicop::create(2, Eigen::VectorXd::Zero(1), 0));
        EXPECT_ANY_THROW(Bicop::create(3, Eigen::VectorXd::Zero(2), 0));
        EXPECT_ANY_THROW(Bicop::create(4, Eigen::VectorXd::Zero(0), 0));
        EXPECT_ANY_THROW(Bicop::create(5, Eigen::VectorXd::Zero(2), 0));
        EXPECT_ANY_THROW(Bicop::create(6, Eigen::VectorXd::Zero(0), 0));
        EXPECT_ANY_THROW(Bicop::create(1001, Eigen::VectorXd::Zero(1), 0));
    }

    TEST(bicop_class, catches_parameters_out_of_bounds) {
        EXPECT_ANY_THROW(Bicop::create(1, Eigen::VectorXd::Constant(1, -1.1), 0));
        EXPECT_ANY_THROW(Bicop::create(2, Eigen::VectorXd::Zero(2), 0));
        EXPECT_ANY_THROW(Bicop::create(2, Eigen::VectorXd::Constant(2, 4.0), 0));
        EXPECT_ANY_THROW(Bicop::create(3, Eigen::VectorXd::Constant(1, -0.1), 0));
        EXPECT_ANY_THROW(Bicop::create(4, Eigen::VectorXd::Constant(1, 1000.0), 0));
        EXPECT_ANY_THROW(Bicop::create(5, Eigen::VectorXd::Constant(1, 10000.0), 0));
        EXPECT_ANY_THROW(Bicop::create(6, Eigen::VectorXd::Constant(1, -0.1), 0));
    }

    TEST(bicop_class, catches_wrong_rotation) {
        EXPECT_ANY_THROW(Bicop::create(1, Eigen::VectorXd::Zero(1), -10));
        EXPECT_ANY_THROW(Bicop::create(1, Eigen::VectorXd::Zero(1), 10));
    }
}