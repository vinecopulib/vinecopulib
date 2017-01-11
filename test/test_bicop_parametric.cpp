/*
    Copyright 2016 Thibault Vatter, Thomas Nagler

    This file is part of vinecoplib.

    vinecoplib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecoplib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
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