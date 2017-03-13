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

#include "src_test/include/parbicop_test.hpp"
#include "src_test/include/test_tools.hpp"
#include "rscript.hpp"

namespace {
    // Test if the C++ implementation of the par_to_tau and tau_to_par is correct
    TYPED_TEST(ParBicopTest, par_to_tau_is_correct) {
        this->setup_parameters();
        std::string command = std::string(RSCRIPT) + "../test/test_bicop_parametric.R";
        command = command + " "  + std::to_string(this->get_n());
        command = command + " "  + std::to_string(this->get_family());
        command = command + " "  + std::to_string(this->get_par());
        command = command + " "  + std::to_string(this->get_par2());
        system(command.c_str());

        if (this->needs_check_) {
            MatXd results = read_matxd("temp");
            VecXd par = this->par_bicop_.get_parameters();
            ASSERT_TRUE(fabs(this->par_bicop_.parameters_to_tau(par) - results(0,0)) < 1e-4);
        }
    }

    // Test if the C++ implementation of the PDF is correct
    TYPED_TEST(ParBicopTest, pdf_is_correct) {
        this->setup_parameters();
        int n = this->get_n();
        if (this->needs_check_) {
            MatXd results = read_matxd("temp");

            // evaluate in C++
            VecXd f = this->par_bicop_.pdf(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,3,n,1), 1e-4));
        }
    }


    // Test if the C++ implementation of the hfunc1 is correct
    TYPED_TEST(ParBicopTest, hfunc1_is_correct) {
        this->setup_parameters();
        int n = this->get_n();
        if (this->needs_check_) {
            MatXd results = read_matxd("temp");

            // evaluate in C++
            VecXd f = this->par_bicop_.hfunc1(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,4,n,1), 1e-4));
        }
    }

    // Test if the C++ implementation of the hfunc2 is correct
    TYPED_TEST(ParBicopTest, hfunc2_is_correct) {
        this->setup_parameters();
        int n = this->get_n();
        if (this->needs_check_) {
            MatXd results = read_matxd("temp");

            // evaluate in C++
            VecXd f = this->par_bicop_.hfunc2(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,5,n,1), 1e-4));
        }
    }

    // Test if the C++ implementation of the hinv1 is correct
    TYPED_TEST(ParBicopTest, hinv1_is_correct) {
        this->setup_parameters();
        int n = this->get_n();
        if (this->needs_check_) {
            MatXd results = read_matxd("temp");

            // evaluate in C++
            VecXd f = this->par_bicop_.hinv1(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,6,n,1), 1e-4));
        }
    }

    // Test if the C++ implementation of the hinv2 is correct
    TYPED_TEST(ParBicopTest, hinv2_is_correct) {
        this->setup_parameters();
        int n = this->get_n();
        if (this->needs_check_) {
            MatXd results = read_matxd("temp");

            // evaluate in C++
            VecXd f = this->par_bicop_.hinv2(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,7,n,1), 1e-4));
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
