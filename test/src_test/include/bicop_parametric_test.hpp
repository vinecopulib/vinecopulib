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

#pragma once

#include "gtest/gtest.h"
#include "include/bicop.hpp"
#include "r_instance.hpp"

// Fake test class allowing access to the R instance
class FakeParBicopTest : public ::testing::Test {
public:
    // all the methods are applied to the R instance
    VecXd eval_in_R(RInstance *rinstance_ptr, std::string eval_fct, int start); // evaluate a function in R
    void set_family(RInstance *rinstance_ptr, int family);
    void set_parameters(RInstance *rinstance_ptr, VecXd parameters);
    void set_rotation(RInstance *rinstance_ptr, int rotation);
    int get_family(RInstance *rinstance_ptr);
    VecXd get_parameters(RInstance *rinstance_ptr);
    int get_rotation(RInstance *rinstance_ptr);
    int get_n(RInstance *rinstance_ptr);
    void change_n(RInstance *rinstance_ptr, int n);
    double get_tau(RInstance *rinstance_ptr);
    void set_tau(RInstance *rinstance_ptr, double tau);
    MatXd get_U(RInstance *rinstance_ptr);
};

// A ParBicop fixture class template
template <typename T>
class ParBicopTest : public FakeParBicopTest {
public:
    void setup_parameters(RInstance *rinstance_ptr) {
        // set family and extracts kendall's tau from the R instance
        int family = this->par_bicop_.get_family();
        this->set_family(rinstance_ptr, family);
        double tau = this->get_tau(rinstance_ptr);

        // set the rotion and compute the parameters vector
        this->par_bicop_.set_rotation(rinstance_ptr->get_rotation());
        VecXd parameters = this->par_bicop_.get_parameters();
        if (parameters.size() < 2)
        {
            parameters = this->par_bicop_.tau_to_parameters(tau);
        }
        else
        {
            if (family == 2)
            {
                parameters = this->par_bicop_.tau_to_parameters(tau);
                parameters(1) = 4;
            }
            else if (family == 7)
            {
                parameters(1) = 1.5;
                parameters(0) = -((2*(1-parameters(1)+parameters(1)*std::fabs(tau)))/(parameters(1)*(-1+std::fabs(tau))));
            }
            else
            {
                double delta = 1.5;
                if (family == 10)
                {
                    delta = 0.8;
                }
                VecXd tau_v = VecXd::Constant(1,std::fabs(tau));
                auto f = [this, delta](const VecXd &v) {
                    VecXd par = VecXd::Constant(2, delta);
                    par(0) = v(0);
                    return VecXd::Constant(1, std::fabs(this->par_bicop_.parameters_to_tau(par)));
                };
                parameters(0) = invert_f(tau_v, f, 1+1e-6, 100)(0);
                parameters(1) = delta;
            }
        }
        // set the parameters vector for the ParBicop and R instance
        this->par_bicop_.set_parameters(parameters);
        this->set_parameters(rinstance_ptr, this->par_bicop_.get_parameters());

        // whether checks need to be done
        needs_check_ = true;
        std::vector<int> rotation_less_fams = {0, 1, 2, 5};
        bool is_rotless = is_member(this->par_bicop_.get_family(), rotation_less_fams);
        if (is_rotless)
            needs_check_ = (rinstance_ptr->get_rotation() == 0);
    }
protected:
    ParBicopTest() : par_bicop_() {}
    T par_bicop_;
    bool needs_check_;
};

// Create a list of types, each of which will be used as the test fixture's 'T'
typedef ::testing::Types<IndepBicop, GaussBicop, StudentBicop, ClaytonBicop, GumbelBicop, FrankBicop, JoeBicop,
        Bb1Bicop, Bb6Bicop, Bb7Bicop, Bb8Bicop> ParBicopTypes;
TYPED_TEST_CASE(ParBicopTest, ParBicopTypes);

