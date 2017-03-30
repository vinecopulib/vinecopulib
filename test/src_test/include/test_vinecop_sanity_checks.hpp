// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "include/vinecop/class.hpp"

namespace test_vinecop_sanity_checks {
    using namespace vinecopulib;

    TEST(vinecop_sanity_checks, catches_wrong_tree) {
        Vinecop vinecop(3);
        EXPECT_ANY_THROW(vinecop.get_pair_copula(3, 0));
    }

    TEST(vinecop_sanity_checks, catches_wrong_edge) {
        Vinecop vinecop(3);
        EXPECT_ANY_THROW(vinecop.get_pair_copula(0, -1));
        EXPECT_ANY_THROW(vinecop.get_pair_copula(0, 2));
    }

    TEST(vinecop_sanity_checks, catches_wrong_size) {
        auto pair_copulas = Vinecop::make_pair_copula_store(3);
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(3,3);
        mat << 2,2,2,1,1,0,3,0,0;
        Vinecop vinecop(pair_copulas, mat);

        Eigen::MatrixXd U = Eigen::MatrixXd::Constant(4,4,0.5);
        EXPECT_ANY_THROW(vinecop.pdf(U));
        EXPECT_ANY_THROW(vinecop.inverse_rosenblatt(U));
        EXPECT_ANY_THROW(vinecop.select_all(U));
        EXPECT_ANY_THROW(vinecop.select_families(U));

        pair_copulas.pop_back();
        EXPECT_ANY_THROW(Vinecop(pair_copulas, mat));
        pair_copulas = Vinecop::make_pair_copula_store(3);
        pair_copulas[0].pop_back();
        EXPECT_ANY_THROW(Vinecop(pair_copulas, mat));
    }
}
