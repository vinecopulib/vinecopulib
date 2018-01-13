// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "kernel_test.hpp"

namespace test_bicop_kernel {
using namespace vinecopulib;

TEST_P(TrafokernelTest, trafo_kernel_sanity_checks) {
    auto values = bicop_.get_parameters();
    EXPECT_ANY_THROW(bicop_.set_parameters(values.block(0, 0, 30, 1)));
    EXPECT_ANY_THROW(bicop_.set_parameters(values.block(0, 0, 1, 30)));
    EXPECT_ANY_THROW(bicop_.set_parameters(-1 * values));
}

TEST_P(TrafokernelTest, trafo_kernel_fit) {
    bicop_.fit(u, controls);
}

TEST_P(TrafokernelTest, trafo_kernel_eval_funcs) {
    bicop_.fit(u, controls);

    EXPECT_GE(bicop_.pdf(u).minCoeff(), 0.0);
    EXPECT_GE(bicop_.cdf(u).minCoeff(), 0.0);
    EXPECT_GE(bicop_.hfunc1(u).minCoeff(), 0.0);
    EXPECT_GE(bicop_.hfunc2(u).minCoeff(), 0.0);
    EXPECT_GE(bicop_.hinv1(u).minCoeff(), 0.0);
    EXPECT_GE(bicop_.hinv2(u).minCoeff(), 0.0);
    EXPECT_LE(bicop_.cdf(u).maxCoeff(), 1.0);
    EXPECT_LE(bicop_.hfunc1(u).maxCoeff(), 1.0);
    EXPECT_LE(bicop_.hfunc2(u).maxCoeff(), 1.0);
    EXPECT_LE(bicop_.hinv1(u).maxCoeff(), 1.0);
    EXPECT_LE(bicop_.hinv2(u).maxCoeff(), 1.0);
    EXPECT_GE(bicop_.calculate_npars(), 0.0);
    EXPECT_LE(bicop_.calculate_npars(), 100.0);

    u(0, 0) = std::numeric_limits<double>::quiet_NaN();
    u(1, 1) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_NO_THROW(bicop_.pdf(u.block(0, 0, 10, 2)));
    EXPECT_TRUE(bicop_.pdf(u.block(0, 0, 1, 2)).array().isNaN()(0));
    EXPECT_NO_THROW(bicop_.cdf(u.block(0, 0, 10, 2)));
    EXPECT_TRUE(bicop_.cdf(u.block(0, 0, 1, 2)).array().isNaN()(0));
    EXPECT_NO_THROW(bicop_.hfunc1(u.block(0, 0, 10, 2)));
    EXPECT_TRUE(bicop_.hfunc1(u.block(0, 0, 1, 2)).array().isNaN()(0));
    EXPECT_NO_THROW(bicop_.hinv1(u.block(0, 0, 10, 2)));
    EXPECT_TRUE(bicop_.hinv1(u.block(0, 0, 1, 2)).array().isNaN()(0));
    EXPECT_NO_THROW(bicop_.hfunc2(u.block(0, 0, 10, 2)));
    EXPECT_TRUE(bicop_.hfunc2(u.block(0, 0, 1, 2)).array().isNaN()(0));
    EXPECT_NO_THROW(bicop_.hinv2(u.block(0, 0, 10, 2)));
    EXPECT_TRUE(bicop_.hinv2(u.block(0, 0, 1, 2)).array().isNaN()(0));
    EXPECT_NO_THROW(bicop_.loglik(u.block(0, 0, 10, 2)));
}

TEST_P(TrafokernelTest, trafo_kernel_select) {
    auto newcop = Bicop(u, controls);
    EXPECT_EQ(newcop.get_family(), BicopFamily::tll);
}

TEST_P(TrafokernelTest, trafo_kernel_flip) {
    auto pdf = bicop_.pdf(u);
    u.col(0).swap(u.col(1));
    bicop_.flip();
    auto pdf_flipped = bicop_.pdf(u);
    EXPECT_TRUE(pdf.isApprox(pdf_flipped, 1e-10));
}

TEST_P(TrafokernelTest, trafo_kernel_tau) {
    double tau = bicop_.parameters_to_tau(bicop_.get_parameters());
    EXPECT_GE(tau, -1.0);
    EXPECT_LE(tau, 1.0);
}


INSTANTIATE_TEST_CASE_P(
    TrafokernelTest,
    TrafokernelTest,
    ::testing::Values("constant", "linear", "quadratic")
);
}
