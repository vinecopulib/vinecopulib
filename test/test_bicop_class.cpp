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

namespace {
    // Test if the C++ implementation of the par_to_tau and tau_to_par is correct
    TEST(bicop_class, creates_right_copula) {
        BicopPtr bicop = Bicop::create(1, VecXd::Ones(1), 90);
        EXPECT_EQ(bicop->get_family(), 1);
        EXPECT_EQ(bicop->get_rotation(), 90);
        EXPECT_EQ(bicop->get_parameters(), VecXd::Ones(1));
    }
    
    TEST(bicop_class, catches_wrong_parameter_size) {
        EXPECT_ANY_THROW(Bicop::create(0, VecXd::Zero(1), 0));
        EXPECT_ANY_THROW(Bicop::create(1, VecXd::Zero(0), 0));
        EXPECT_ANY_THROW(Bicop::create(2, VecXd::Zero(1), 0));
        EXPECT_ANY_THROW(Bicop::create(3, VecXd::Zero(2), 0));
        EXPECT_ANY_THROW(Bicop::create(4, VecXd::Zero(0), 0));
        EXPECT_ANY_THROW(Bicop::create(5, VecXd::Zero(2), 0));
        EXPECT_ANY_THROW(Bicop::create(6, VecXd::Zero(0), 0));
        EXPECT_ANY_THROW(Bicop::create(1001, VecXd::Zero(1), 0));
    }
    
    TEST(bicop_class, catches_parameters_out_of_bounds) {
        EXPECT_ANY_THROW(Bicop::create(1, VecXd::Constant(1, -1.1), 0));
        EXPECT_ANY_THROW(Bicop::create(2, VecXd::Zero(2), 0));
        EXPECT_ANY_THROW(Bicop::create(2, VecXd::Constant(2, 4.0), 0));
        EXPECT_ANY_THROW(Bicop::create(3, VecXd::Constant(1, -0.1), 0));
        EXPECT_ANY_THROW(Bicop::create(4, VecXd::Constant(1, 1000.0), 0));
        EXPECT_ANY_THROW(Bicop::create(5, VecXd::Constant(1, 10000.0), 0));
        EXPECT_ANY_THROW(Bicop::create(6, VecXd::Constant(1, -0.1), 0));
    }
    
    TEST(bicop_class, catches_wrong_rotation) {
        EXPECT_ANY_THROW(Bicop::create(1, VecXd::Zero(1), -10));
        EXPECT_ANY_THROW(Bicop::create(1, VecXd::Zero(1), 10));
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
