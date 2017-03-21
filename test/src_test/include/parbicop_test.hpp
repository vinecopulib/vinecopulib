// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "include/bicop.hpp"
#include "include/tools_stl.hpp"
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

using namespace vinecopulib;

// Fake test class allowing access to the R instance
class FakeParBicopTest : public ::testing::Test {
public:
    void set_family(BicopFamily family, int rotation);
    void set_parameters(Eigen::VectorXd parameters);
    void set_n(int n);
    int get_n();
    int get_family();
    double get_par();
    double get_par2();
private:
    int n_;
    int family_;
    double par_;
    double par2_;
};

// A ParBicop fixture class template
template <typename T>
class ParBicopTest : public FakeParBicopTest {
public:
    void setup_parameters(int rotation = 0) {
        int n = (int) 5e3;
        double tau = 0.5; // should be positive
        auto family = par_bicop_.get_family();
        bicop_ = Bicop::create(family, rotation);

        this->set_family(family, rotation);
        this->set_n(n);

        auto parameters = bicop_->get_parameters();
        if (parameters.size() < 2) {
            parameters = bicop_->tau_to_parameters(tau); 
        } else {
            if (family == BicopFamily::student) {
                parameters(1) = 4;
            } else if (family == BicopFamily::bb1) {
                parameters(1) = 1.5;
                parameters(0) = -((2*(1-parameters(1)+parameters(1)*tau))/(parameters(1)*(-1+tau)));
            } else {
                double delta = 1.5;
                if (family == BicopFamily::bb8)
                    delta = 0.8;
                auto tau_v = Eigen::VectorXd::Constant(1, std::fabs(tau));
                auto f = [this, delta](const Eigen::VectorXd &v) {
                    Eigen::VectorXd par(2);
                    par << v(0), delta;
                    auto tau = bicop_->parameters_to_tau(par);
                    return Eigen::VectorXd::Constant(1, std::fabs(tau));
                };
                parameters(0) = invert_f(tau_v, f, 1 + 1e-6, 100)(0);
                parameters(1) = delta;
            }
        }
        // set the parameters vector for the ParBicop
        bicop_->set_parameters(parameters);

        // whether checks need to be done and deal with the rotation for VineCopula
        needs_check_ = true;
        if (tools_stl::is_member(family, bicop_families::rotationless)) {
            needs_check_ = (rotation == 0);
        } else {
            if (tools_stl::is_member(rotation, {90, 270}))
                parameters *= -1;
        }

        // set the parameters vector for R
        this->set_parameters(parameters);
    }
    
protected:
    ParBicopTest() : par_bicop_() {}
    T par_bicop_;
    std::shared_ptr<Bicop> bicop_;
    bool needs_check_;
};

// Create a list of types, each of which will be used as the test fixture's 'T'
typedef ::testing::Types<
    IndepBicop, GaussianBicop, StudentBicop, ClaytonBicop, GumbelBicop,
    FrankBicop, JoeBicop, Bb1Bicop, Bb6Bicop, Bb7Bicop, Bb8Bicop
> ParBicopTypes;
TYPED_TEST_CASE(ParBicopTest, ParBicopTypes);
