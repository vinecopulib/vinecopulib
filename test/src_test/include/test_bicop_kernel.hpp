// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "kernel_test.hpp"

namespace test_bicop_trafokernel {
    using namespace vinecopulib;

    TEST_F(TrafokernelTest, trafo_kernel_fit) {
        bicop_.fit(u, std::string(""));
    }

    TEST_F(TrafokernelTest, trafo_kernel_eval_funcs) {
        bicop_.fit(u, std::string(""));

        EXPECT_GE(bicop_.pdf(u).minCoeff(), 0.0);

        EXPECT_GE(bicop_.hfunc1(u).minCoeff(), 0.0);
        EXPECT_GE(bicop_.hfunc2(u).minCoeff(), 0.0);
        EXPECT_GE(bicop_.hinv1(u).minCoeff(), 0.0);
        EXPECT_GE(bicop_.hinv2(u).minCoeff(), 0.0);

        EXPECT_LE(bicop_.hfunc1(u).maxCoeff(), 1.0);
        EXPECT_LE(bicop_.hfunc2(u).maxCoeff(), 1.0);
        EXPECT_LE(bicop_.hinv1(u).maxCoeff(), 1.0);
        EXPECT_LE(bicop_.hinv2(u).maxCoeff(), 1.0);

        EXPECT_GE(bicop_.calculate_npars(), 0.0);
        EXPECT_LE(bicop_.calculate_npars(), 100.0);
    }

    TEST_F(TrafokernelTest, trafo_kernel_select) {
        auto newcop = Bicop(u, {BicopFamily::tll0});
        EXPECT_EQ(newcop.get_family(), BicopFamily::tll0);
    }

    TEST_F(TrafokernelTest, trafo_kernel_flip) {
        auto pdf = bicop_.pdf(u);
        u.col(0).swap(u.col(1));
        bicop_.flip();
        auto pdf_flipped = bicop_.pdf(u);
        EXPECT_TRUE(pdf.isApprox(pdf_flipped, 1e-10));
    }
}
