// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>


namespace test_bicop_discrete {

using namespace vinecopulib;

TEST(bicop_discrete, rotations_are_correct) {
    for (auto rot : {0, 90, 180, 270}) {
        auto bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
        auto u = bc.simulate(1000, true, {1});

        Eigen::MatrixXd u_new(u.rows(), 4);
        u_new.block(0, 0, u.rows(), 2) = u;
        u_new.block(0, 2, u.rows(), 2) = u;

        u_new.col(0) = (u.col(0).array() * 2).ceil() / 2;
        u_new.col(2) = (u.col(0).array() * 2).floor() / 2;
        bc.set_discrete_vars({1});
        EXPECT_GE(bc.pdf(u_new).minCoeff(), 0);
        bc.fit(u_new);
        EXPECT_NEAR(bc.get_parameters()(0), 3, 0.3);

        u_new.block(0, 0, u.rows(), 2) = u;
        u_new.block(0, 2, u.rows(), 2) = u;
        u_new.col(1) = (u.col(1).array() * 2).ceil() / 2;
        u_new.col(3) = (u.col(1).array() * 2).floor() / 2;
        bc.set_discrete_vars({2});
        EXPECT_GE(bc.pdf(u_new).minCoeff(), 0);
        bc.fit(u_new);
        EXPECT_NEAR(bc.get_parameters()(0), 3, 0.3);

        u_new.col(0) = (u_new.col(0).array() * 2).ceil() / 2.0;
        u_new.col(2) = (u_new.col(2).array() * 2).floor() / 2.0;
        bc.set_discrete_vars({1, 2});
        EXPECT_GE(bc.pdf(u_new).minCoeff(), 0);
        bc.fit(u_new);
        EXPECT_NEAR(bc.get_parameters()(0), 3, 0.3);
    }
}

TEST(bicop_discrete, playground) {
    for (auto rot : {0}) {
        auto bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
        auto u = bc.simulate(20, true, {1});
        // std::cout << u << std::endl<< std::endl;

        Eigen::MatrixXd u_new(u.rows(), 4);
        u_new.block(0, 0, u.rows(), 2) = u;
        u_new.block(0, 2, u.rows(), 2) = u;

        u_new.col(0) = (u.col(0).array() * 2).ceil() / 2;
        u_new.col(2) = (u.col(0).array() * 2).floor() / 2;
        u_new.col(1) = (u.col(1).array() * 2).ceil() / 2;
        u_new.col(3) = (u.col(1).array() * 2).floor() / 2;
        bc.set_discrete_vars({1, 2});
        bc.fit(u_new);
    }
}

}
