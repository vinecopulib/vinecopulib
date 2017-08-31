// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_integration.hpp>
#include <cmath>

using namespace vinecopulib;
namespace test_tools_integration {

    TEST(test_tools_stats, legendre_rule_is_correct) {

        // Get quadrature nodes and weights
        size_t N = 25;
        Eigen::MatrixXd nodes_weights = tools_integration::legendre_rule(N);
        Eigen::VectorXd nodes = nodes_weights.col(0);
        Eigen::VectorXd weights = nodes_weights.col(1);

        // Mean and variance of uniform distribution
        double m1 = nodes.cwiseProduct(weights).sum();
        EXPECT_NEAR(m1, 0.5, 1e-4);
        double m2 = nodes.cwiseAbs2().cwiseProduct(weights).sum();
        EXPECT_NEAR(m2 - m1 * m1, 1/12.0, 1e-4);

        // int^1_0 exp(x) dx
        auto f = [](const double x, const double w) { return w * std::exp(x); };
        double m3 = nodes.binaryExpr(weights, f).sum();
        EXPECT_NEAR(m3, std::exp(1.0) - 1.0, 1e-4);

    }

    TEST(test_tools_stats, integrate_zero_to_one_is_correct) {

        // int^1_0 exp(x) dx
        auto f = [](const double x) { return std::exp(x); };
        double m1 = tools_integration::integrate_zero_to_one(f);
        EXPECT_NEAR(m1, std::exp(1.0) - 1.0, 1e-4);
    }
}
