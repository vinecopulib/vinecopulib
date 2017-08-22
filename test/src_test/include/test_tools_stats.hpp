// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

using namespace vinecopulib;

namespace test_tools_stats {

    TEST(test_tools_stats, to_pseudo_obs_is_correct) {

        int n = 9;

        // X1 = (1,...,n) and X2 = (n, ..., 1)
        // X = (X1, X2)
        Eigen::MatrixXd X(n,2);
        X.col(0) = Eigen::VectorXd::LinSpaced(n,1,n);
        X.col(1) = Eigen::VectorXd::LinSpaced(n,n,1);

        // U = pobs(X)
        Eigen::MatrixXd U = tools_stats::to_pseudo_obs(X);
        for (int i = 0; i < 9; i++) {
            EXPECT_NEAR(U(i,0), (i+1.0)*0.1, 1e-2);
            EXPECT_NEAR(U(i,1), 1.0-(i+1.0)*0.1, 1e-2);
        }
    }

    TEST(test_tools_stats, pairwise_dep_measures_are_correct) {    
        // Independence
        Eigen::Matrix<double, Eigen::Dynamic, 2> X(10, 2);
        X <<  4,  6,
              7,  3,
              1,  8,
              9, 10,
             10,  5,
              6,  7,
              5,  9,
              8,  1,
              3,  2,
              2,  4;
        X = X.array() / 10.0;
        EXPECT_NEAR(tools_stats::pairwise_tau(X),    0.0, 5e-1);
        EXPECT_NEAR(tools_stats::pairwise_cor(X),    0.0, 5e-1);
        EXPECT_NEAR(tools_stats::pairwise_rho(X),    0.0, 5e-1);
        EXPECT_NEAR(tools_stats::pairwise_hoeffd(X), 0.0, 5e-1);
        
        // Perfect (negative) dependence        
        X.col(1) = 1.0 - X.col(0).array();
        EXPECT_NEAR(tools_stats::pairwise_tau(X),   -1.0, 1e-3);
        EXPECT_NEAR(tools_stats::pairwise_cor(X),   -1.0, 1e-3);
        EXPECT_NEAR(tools_stats::pairwise_rho(X),   -1.0, 1e-3);
        EXPECT_NEAR(tools_stats::pairwise_hoeffd(X), 1.0, 1e-3);


    }
    
    TEST(test_tools_stats, dependence_matrix_works) {
        auto u = tools_stats::simulate_uniform(10, 4);
        auto mat = tools_stats::dependence_matrix(u, "tau");
        EXPECT_TRUE((mat.diagonal().array() == 1.0).all());
        EXPECT_TRUE(mat(0, 1) == mat(1, 0));

        // only check if it works from here one
        tools_stats::dependence_matrix(u, "cor");
        tools_stats::dependence_matrix(u, "rho");
        tools_stats::dependence_matrix(u, "hoeffd");
        EXPECT_ANY_THROW(tools_stats::dependence_matrix(u, "other"));
    }

    TEST(test_tools_stats, ghalton_is_correct) {

        size_t d = 2;
        size_t n = 1e1;
        size_t N = 1e4;

        auto cop = Bicop(BicopFamily::gaussian);
        auto u = cop.simulate(n);
        auto U = tools_stats::ghalton(N, d);
        auto U2 = tools_stats::simulate_uniform(N, d);

        Eigen::VectorXd x(N), p(n), x2(N), p2(n);
        p2 = Eigen::VectorXd::Zero(n);
        for (size_t i = 0; i < n; i++) {
            auto f = [i, u](const double& u1, const double& u2) {
                return (u1 <= u(i,0) && u2 <= u(i,1)) ? 1.0 : 0.0;
            };
            x = U.col(0).binaryExpr(cop.hinv1(U), f);
            p(i) = x.sum()/N;
            x2 = U2.col(0).binaryExpr(cop.hinv1(U2), f);
            p2(i) = x2.sum()/N;
        }

        if (p2.isApprox(cop.cdf(u), 1e-2))
            ASSERT_TRUE(p.isApprox(cop.cdf(u), 1e-2));

    }
}
