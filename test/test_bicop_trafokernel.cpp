/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
