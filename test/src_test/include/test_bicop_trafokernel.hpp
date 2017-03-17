// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "trafokernel_test.hpp"

namespace test_bicop_trafokernel {
    using namespace vinecopulib;

    TEST_F(TrafokernelTest, trafo_kernel_fit) {
        cop->fit(u, std::string(""));
    }

    TEST_F(TrafokernelTest, trafo_kernel_eval_funcs) {
        cop->fit(u, std::string(""));

        EXPECT_GE(cop->pdf(u).minCoeff(), 0.0);

        EXPECT_GE(cop->hfunc1(u).minCoeff(), 0.0);
        EXPECT_GE(cop->hfunc2(u).minCoeff(), 0.0);
        EXPECT_GE(cop->hinv1(u).minCoeff(), 0.0);
        EXPECT_GE(cop->hinv2(u).minCoeff(), 0.0);

        EXPECT_LE(cop->hfunc1(u).maxCoeff(), 1.0);
        EXPECT_LE(cop->hfunc2(u).maxCoeff(), 1.0);
        EXPECT_LE(cop->hinv1(u).maxCoeff(), 1.0);
        EXPECT_LE(cop->hinv2(u).maxCoeff(), 1.0);

        EXPECT_GE(cop->calculate_npars(), 0.0);
        EXPECT_LE(cop->calculate_npars(), 100.0);
    }

    TEST_F(TrafokernelTest, trafo_kernel_select) {
        auto newcop = Bicop::select(u, {BicopFamily::tll0});
        EXPECT_EQ(newcop->get_family(), BicopFamily::tll0);
    }

    TEST_F(TrafokernelTest, trafo_kernel_flip) {
        auto pdf = cop->pdf(u);
        u.col(0).swap(u.col(1));
        cop->flip();
        auto pdf_flipped = cop->pdf(u);
        EXPECT_TRUE(pdf.isApprox(pdf_flipped, 1e-10));
    }
}
