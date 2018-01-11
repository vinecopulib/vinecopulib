// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>

namespace test_bicop_select {
using namespace vinecopulib;

TEST(bicop_select, works_in_parallel) {
    Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
    auto u = cop.simulate(15);
    Bicop fit1, fit2;
    fit1.select(u);
    FitControlsBicop controls;
    controls.set_num_threads(2);
    fit2.select(u, controls);
    EXPECT_EQ(fit1.get_family(), fit2.get_family());
    EXPECT_EQ(fit1.get_parameters(), fit2.get_parameters());
}

TEST(bicop_select, allows_all_selcrits) {
    Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
    auto u = cop.simulate(15);
    FitControlsBicop controls;
    controls.set_selection_criterion("loglik");
    cop.select(u, controls);
    controls.set_selection_criterion("aic");
    cop.select(u, controls);
    controls.set_selection_criterion("bic");
    cop.select(u, controls);
}

}
