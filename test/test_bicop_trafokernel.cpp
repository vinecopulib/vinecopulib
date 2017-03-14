/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

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
        auto cop = Bicop::select(u, "bic", {1001}, true, true, "mle");
        EXPECT_EQ(cop->get_family(), 1001);
    }

    TEST(KernelBicoptest, trafo_kernel_flip) {
        auto u = simulate_uniform(100, 2);
        auto cop = Bicop::select(u, "bic", {1001}, true, true, "mle");
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
