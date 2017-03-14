// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include "include/bicop.hpp"
#include <iostream>

namespace {

    TEST(KernelBicopTest, trafo_kernel_fit) {
        auto cop = Bicop::create(1001, 0);
        auto u = simulate_uniform(100, 2);
        cop->fit(u, std::string(""));
    }

    TEST(KernelBicop, trafo_kernel_eval_funcs) {
        auto cop = Bicop::create(1001, 0);
        auto u = simulate_uniform(100, 2);
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

    TEST(KernelBicoptest, trafo_kernel_select) {
        auto u = simulate_uniform(100, 2);
        auto cop = Bicop::select(u, {1001});
        EXPECT_EQ(cop->get_family(), 1001);
    }

    TEST(KernelBicoptest, trafo_kernel_flip) {
        auto u = simulate_uniform(100, 2);
        auto cop = Bicop::select(u, {1001});
        auto pdf = cop->pdf(u);
        u.col(0).swap(u.col(1));
        cop->flip();
        auto pdf_flipped = cop->pdf(u);
        EXPECT_TRUE(pdf.isApprox(pdf_flipped, 1e-10));
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
