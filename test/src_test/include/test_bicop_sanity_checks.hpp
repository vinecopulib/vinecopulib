// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>

namespace test_bicop_sanity_checks {
using namespace vinecopulib;

TEST(bicop_sanity_checks, catches_wrong_parameter_size) {
    for (auto family : bicop_families::nonparametric) {
        EXPECT_ANY_THROW(Bicop(family, 0, Eigen::VectorXd::Zero(1)));
    }
    for (auto family : bicop_families::one_par) {
        EXPECT_ANY_THROW(Bicop(family, 0, Eigen::VectorXd::Zero(2)));
    }
    for (auto family : bicop_families::two_par) {
        EXPECT_ANY_THROW(Bicop(family, 0, Eigen::VectorXd::Zero(1)));
    }
}

TEST(bicop_sanity_checks, catches_parameters_out_of_bounds) {
    auto cop = Bicop(BicopFamily::gaussian);
    auto wrong_par = Eigen::VectorXd::Constant(1, 1.01);
    EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 0, wrong_par));
    EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 0, -wrong_par));
}

TEST(bicop_sanity_checks, catches_wrong_rotation) {
    EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, -10));
    EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 10));
}

TEST(bicop_sanity_checks, select_can_handle_zeros_and_ones) {
    Bicop bicop;
    auto u = bicop.simulate(10);
    u(0, 0) = 0.0;
    u(1, 0) = 1.0;
    EXPECT_NO_THROW(bicop.select(u));
}

}
