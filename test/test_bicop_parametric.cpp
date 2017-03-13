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

#include "src_test/include/bicop_parametric_test.hpp"

RInstance *rinstance_ptr = new RInstance;

namespace {
    // Test if the C++ implementation of the par_to_tau and tau_to_par is correct
    TYPED_TEST(ParBicopTest, par_to_tau_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "BiCopPar2Tau(,,)";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 13);

            // assert approximate equality
            VecXd par = this->get_parameters(rinstance_ptr);
            ASSERT_TRUE(fabs(this->par_bicop_.parameters_to_tau(par) - f(0)) < 1e-4);
        }
    }

    // Test if the C++ implementation of the PDF is correct
    TYPED_TEST(ParBicopTest, pdf_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "BiCopPDF(u1,u2,,,)";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 15);

            // evaluate in C++
            VecXd f2 = this->par_bicop_.pdf(U);

            // assert approximate equality
            MatXd output(f2.size(), 2);
            output.col(0) = f;
            output.col(1) = f2;
            ASSERT_TRUE(f.isApprox(f2, 1e-4));
        }
    }

    // Test if the C++ implementation of the hfunc1 is correct
    TYPED_TEST(ParBicopTest, hfunc1_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "BiCopHfunc1(u1,u2,,,)";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 18);

            // evaluate in C++
            VecXd f2 = this->par_bicop_.hfunc1(U);

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(f2, 1e-4)) << f(0) << " != " << f2(0);
        }
    }

    // Test if the C++ implementation of the hfunc2 is correct
    TYPED_TEST(ParBicopTest, hfunc2_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "BiCopHfunc2(u1,u2,,,)";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 18);

            // evaluate in C++
            VecXd f2 = this->par_bicop_.hfunc2(U);

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(f2, 1e-4)) << f(0) << " != " << f2(0);
        }
    }

    // Test if the C++ implementation of the hinv1 is correct
    TYPED_TEST(ParBicopTest, hinv1_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "BiCopHinv1(u1,u2,,,)";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 17);

            // evaluate in C++
            VecXd f2 = this->par_bicop_.hinv1(U);

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(f2, 1e-4));
        }
    }

    // Test if the C++ implementation of the hinv2 is correct
    TYPED_TEST(ParBicopTest, hinv2_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "BiCopHinv2(u1,u2,,,)";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 17);

            // evaluate in C++
            VecXd f2 = this->par_bicop_.hinv2(U);

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(f2, 1e-4));
        }
    }

    // Test if the C++ implementation of the loglik is correct
    TYPED_TEST(ParBicopTest, loglik_is_correct) {
        this->setup_parameters(rinstance_ptr);
        if (this->needs_check_) {
            // evaluate in R
            MatXd U = this->get_U(rinstance_ptr);
            std::string eval_fct = "sum(log(BiCopPDF(u1,u2,,,)))";
            VecXd f = this->eval_in_R(rinstance_ptr, eval_fct, 23);

            // evaluate in C++
            double f2 = this->par_bicop_.loglik(U);

            // assert approximate equality
            ASSERT_TRUE(fabs(f2 - f(0)) < 1e-2);
        }
    }
}

int main(int argc, char **argv) {
    rinstance_ptr->set_rotation(0);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
