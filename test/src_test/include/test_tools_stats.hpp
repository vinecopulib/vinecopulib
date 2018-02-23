// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "test_vinecop_sanity_checks.hpp"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace test_tools_stats {

using namespace vinecopulib;

TEST(test_tools_stats, to_pseudo_obs_is_correct) {

    int n = 9;

    // X1 = (1,...,n) and X2 = (n, ..., 1)
    // X = (X1, X2)
    Eigen::MatrixXd X(n, 2);
    X.col(0) = Eigen::VectorXd::LinSpaced(n, 1, n);
    X.col(1) = Eigen::VectorXd::LinSpaced(n, n, 1);

    // U = pobs(X)
    Eigen::MatrixXd U = tools_stats::to_pseudo_obs(X);
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(U(i, 0), (i + 1.0) * 0.1, 1e-2);
        EXPECT_NEAR(U(i, 1), 1.0 - (i + 1.0) * 0.1, 1e-2);
    }

    Eigen::MatrixXd X2 = tools_stats::simulate_uniform(100, 2);
    EXPECT_NO_THROW(tools_stats::to_pseudo_obs(X2, "random"));
    EXPECT_NO_THROW(tools_stats::to_pseudo_obs(X2, "first"));
    EXPECT_ANY_THROW(tools_stats::to_pseudo_obs(X2, "something"));
}

TEST(test_tools_stats, pairwise_dep_measures_are_correct) {
    // Independence
    Eigen::Matrix<double, Eigen::Dynamic, 2> X(10, 2);
    X << 4, 6,
        7, 3,
        1, 8,
        9, 10,
        10, 5,
        6, 7,
        5, 9,
        8, 1,
        3, 2,
        2, 4;
    X = X.array() / 10.0;
    EXPECT_NEAR(tools_stats::pairwise_tau(X), 0.0, 5e-1);
    EXPECT_NEAR(tools_stats::pairwise_cor(X), 0.0, 5e-1);
    EXPECT_NEAR(tools_stats::pairwise_rho(X), 0.0, 5e-1);
    EXPECT_NEAR(tools_stats::pairwise_hoeffd(X), 0.0, 5e-1);

    // Perfect (negative) dependence
    X.col(1) = 1.0 - X.col(0).array();
    EXPECT_NEAR(tools_stats::pairwise_tau(X), -1.0, 1e-3);
    EXPECT_NEAR(tools_stats::pairwise_cor(X), -1.0, 1e-3);
    EXPECT_NEAR(tools_stats::pairwise_rho(X), -1.0, 1e-3);
    EXPECT_NEAR(tools_stats::pairwise_hoeffd(X), 1.0, 1e-3);


}

TEST(test_tools_stats, dependence_matrix_works) {
    auto u = tools_stats::simulate_uniform(10, 4);
    auto mat = tools_stats::dependence_matrix(u, "tau");
    EXPECT_TRUE((mat.diagonal().array() == 1.0).all());
    EXPECT_TRUE(mat(0, 1) == mat(1, 0));

    // only check if it works from here one
    tools_stats::dependence_matrix(u, "cor");
    tools_stats::dependence_matrix(u, "mcor");
    tools_stats::dependence_matrix(u, "rho");
    tools_stats::dependence_matrix(u, "hoeffd");
    EXPECT_ANY_THROW(tools_stats::dependence_matrix(u, "other"));
}

TEST(test_tools_stats, ghalton_is_correct) {

    size_t d = 2;
    size_t n = 10;
    size_t N = 1000;

    auto cop = Bicop(BicopFamily::gaussian);
    auto u = cop.simulate(n);
    auto U = tools_stats::ghalton(N, d);
    auto U2 = tools_stats::simulate_uniform(N, d);

    Eigen::VectorXd x(N), p(n), x2(N), p2(n);
    p2 = Eigen::VectorXd::Zero(n);
    for (size_t i = 0; i < n; i++) {
        auto f = [i, u](const double &u1, const double &u2) {
            return (u1 <= u(i, 0) && u2 <= u(i, 1)) ? 1.0 : 0.0;
        };
        x = U.col(0).binaryExpr(cop.hinv1(U), f);
        p(i) = x.sum() / N;
        x2 = U2.col(0).binaryExpr(cop.hinv1(U2), f);
        p2(i) = x2.sum() / N;
    }

    if (p2.isApprox(cop.cdf(u), 1e-2)) {
        ASSERT_TRUE(p.isApprox(cop.cdf(u), 1e-2));
    }


}

TEST(test_tools_stats, dpqnorm_are_nan_safe) {
    Eigen::VectorXd X = Eigen::VectorXd::Random(10);
    X(0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_NO_THROW(tools_stats::dnorm(X));
    EXPECT_NO_THROW(tools_stats::pnorm(X));
    EXPECT_NO_THROW(tools_stats::qnorm(tools_stats::pnorm(X)));
}

TEST(test_tools_stats, dpt_are_nan_safe) {
    Eigen::VectorXd X = Eigen::VectorXd::Random(10);
    X(0) = std::numeric_limits<double>::quiet_NaN();
    double nu = 4.0;
    EXPECT_NO_THROW(tools_stats::dt(X, nu));
    EXPECT_NO_THROW(tools_stats::pt(X, nu));
    EXPECT_NO_THROW(tools_stats::qt(tools_stats::pt(X, nu), nu));
}

TEST(test_tools_stats, pbvt_and_pbvnorm_are_nan_safe) {
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(10, 2);
    X(0) = std::numeric_limits<double>::quiet_NaN();
    double rho = -0.95;
    int nu = 5;
    EXPECT_NO_THROW(tools_stats::pbvt(X, nu, rho));
    EXPECT_NO_THROW(tools_stats::pbvnorm(X, rho));
}

}
