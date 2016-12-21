/*
    Copyright 2016 Thibault Vatter

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

#ifndef VINECOPLIB_TEST_PAR_BICOP_H
#define VINECOPLIB_TEST_PAR_BICOP_H

#include "gtest/gtest.h"
#include "src/common/include/bicop_indep.hpp"
#include "src/common/include/bicop_normal.hpp"
#include "src/common/include/bicop_student.hpp"
#include "src/common/include/bicop_clayton.hpp"
#include "src/common/include/bicop_gumbel.hpp"
#include "src/common/include/bicop_frank.hpp"
#include "src/common/include/bicop_joe.hpp"
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
        double tau = this->get_tau(rinstance_ptr);
        double par = this->par_bicop_.tau_to_par(tau);
        VecXd parameters = par * VecXd::Ones(2);
        parameters(1) = 4;
        this->par_bicop_.set_parameters(parameters);
        this->set_family(rinstance_ptr, this->par_bicop_.get_family());
        this->set_parameters(rinstance_ptr, this->par_bicop_.get_parameters());
    };
protected:
    ParBicopTest() : par_bicop_() {}
    T par_bicop_;
};

// Create a list of types, each of which will be used as the test fixture's 'T'
typedef ::testing::Types<IndepBicop, NormalBicop, StudentBicop, ClaytonBicop, GumbelBicop, FrankBicop, JoeBicop> ParBicopTypes;
TYPED_TEST_CASE(ParBicopTest, ParBicopTypes);
//typedef ::testing::Types<ClaytonBicop, GumbelBicop> RotatedBicopTypes;
//TYPED_TEST_CASE(ParBicopTest, RotatedBicopTypes);

#endif //VINECOPLIB_TEST_PAR_BICOP_H
