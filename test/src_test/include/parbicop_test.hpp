// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <boost/math/constants/constants.hpp>

using namespace vinecopulib;

// Test class for parametric bivariate copulas
class ParBicopTest
  : public ::testing::TestWithParam<::testing::tuple<BicopFamily, int>>
{
public:
  void set_family(BicopFamily family, int rotation);

  void set_parameters(Eigen::VectorXd parameters);

  int get_n();

  int get_family();

  double get_par();

  double get_par2();

protected:
  int n_;
  int family_;
  double par_;
  double par2_;
  Bicop bicop_;
  bool needs_check_;

  virtual void SetUp()
  {
    n_ = static_cast<int>(5e3);
    auto family = ::testing::get<0>(GetParam());
    auto rotation = ::testing::get<1>(GetParam());
    if (tools_stl::is_member(family, bicop_families::rotationless)) {
      bicop_ = Bicop(family);
    } else {
      bicop_ = Bicop(family, rotation);
    }

    set_family(family, rotation);
    double tau = 0.5; // should be positive
    auto parameters = bicop_.get_parameters();
    if (parameters.size() < 2) {
      parameters = bicop_.tau_to_parameters(tau);
    } else if (parameters.size() == 2) {
      if (family == BicopFamily::student) {
        parameters(0) = sin(tau * boost::math::constants::pi<double>() / 2);
        parameters(1) = 4;
      } else if (family == BicopFamily::bb1) {
        parameters(1) = 1.5;
        parameters(0) = -(2 * (1 - parameters(1) + parameters(1) * tau));
        parameters(0) /= (parameters(1) * (-1 + tau));
      } else {
        double delta = 1.5;
        if (family == BicopFamily::bb8)
          delta = 0.8;
        auto tau_v = Eigen::VectorXd::Constant(1, std::fabs(tau));
        auto f = [this, delta](const Eigen::VectorXd& v) {
          Eigen::VectorXd par(2);
          par(0) = v(0);
          par(1) = delta;
          auto tt = bicop_.parameters_to_tau(par);
          return Eigen::VectorXd::Constant(1, std::fabs(tt));
        };
        parameters(0) = tools_eigen::invert_f(tau_v, f, 1 + 1e-6, 100)(0);
        parameters(1) = delta;
      }
    } else {
      parameters(0) = 0.5;
      parameters(1) = 1;
      parameters(2) = 6;
    }
    // set the parameters vector for the ParBicop
    bicop_.set_parameters(parameters);

    // whether checks need to be done and deal with the rotation for VineCopula
    needs_check_ = true;
    if (tools_stl::is_member(family, bicop_families::rotationless)) {
      needs_check_ = (rotation == 0);
    } else {
      if (tools_stl::is_member(rotation, { 90, 270 })) {
        parameters *= -1;
        if (family == BicopFamily::tawn)
          parameters(0) *= -1;
      }
    }

    // set the parameters vector for R
    set_parameters(parameters);
  }
};
