// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/abstract.hpp>

namespace vinecopulib {

//! @brief An abstract class for parametric copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class ParBicop : public AbstractBicop
{
protected:
  // Getters and setters
  Matrix get_parameters() const;

  Matrix get_parameters_lower_bounds() const;

  Matrix get_parameters_upper_bounds() const;

  void set_parameters(const Matrix& parameters);

  void flip();

  // Data members
  Matrix parameters_;
  Matrix parameters_lower_bounds_;
  Matrix parameters_upper_bounds_;

  void fit(const Matrix& data,
           std::string method,
           double,
           const Vector& weights);

  double get_npars() const;

  void set_npars(const double& npars);

  virtual Vector get_start_parameters(const double tau) = 0;

private:
  double winsorize_tau(double tau) const;

  void adjust_parameters_bounds(Matrix& lb,
                                Matrix& ub,
                                const double& tau,
                                const std::string& method);

  void check_parameters(const Matrix& parameters);

  void check_parameters_size(const Matrix& parameters);

  void check_parameters_upper(const Matrix& parameters);

  void check_parameters_lower(const Matrix& parameters);

  void check_fit_method(const std::string& method);
};
}

#include <vinecopulib/bicop/implementation/parametric.ipp>
