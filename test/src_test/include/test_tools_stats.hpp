// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace test_tools_stats {

    TEST(test_tools_stats, to_pseudo_obs_is_correct) {

        int n = 9;

        // X1 = (1,...,n) and X2 = (n, ..., 1)
        // X = (X1, X2)
        Eigen::MatrixXd X(n,2);
        X.col(0) = Eigen::VectorXd::LinSpaced(n,1,n);
        X.col(1) = Eigen::VectorXd::LinSpaced(n,n,1);

        // U = pobs(X)
        Eigen::MatrixXd U = vinecopulib::tools_stats::to_pseudo_obs(X);
        for (int i = 0; i < 9; i++) {
            EXPECT_NEAR(U(i,0), (i+1.0)*0.1, 1e-2);
            EXPECT_NEAR(U(i,1), 1.0-(i+1.0)*0.1, 1e-2);
        }
    }

    TEST(test_tools_stats, hoeffd_is_correct) {

        int n = 9;

        // X1 = (1,...,n) and X2 = (n, ..., 1)
        // X = (X1, X2)
        Eigen::MatrixXd X(n,2);
        X.col(0) = Eigen::VectorXd::LinSpaced(n,1,n);
        X.col(1) = Eigen::VectorXd::LinSpaced(n,n,1);

        // Perfect dependence
        double computed_hoeffd = vinecopulib::tools_stats::pairwise_hoeffd(X);
        double true_hoeffd = 1.0;
        EXPECT_NEAR(computed_hoeffd, true_hoeffd, 1e-3);

        // Independence
        Eigen::MatrixXd X2(5,2);
        X2 << -2,4,-1,1,0,0,1,2,2,3;
        computed_hoeffd = vinecopulib::tools_stats::pairwise_hoeffd(X2);
        true_hoeffd = 0.0;
        EXPECT_NEAR(computed_hoeffd, true_hoeffd, 1e-3);
    }
}
