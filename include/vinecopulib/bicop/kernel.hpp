// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/abstract.hpp>

namespace vinecopulib {

namespace tools_interpolation {
class InterpolationGrid;
}

//! @brief An abstract class for kernel copulas.
//!
//! Evaluation functions of kernel estimators are implemented efficiently
//! using spline interpolation, see Nagler (2016).
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Nagler, Thomas. *kdecopula: An R Package for the Kernel Estimation of
//! Copula Densities*. arXiv:1603.04229 [stat.CO], 2016
class KernelBicop : public AbstractBicop
{
public:
  KernelBicop();

protected:
  Eigen::VectorXd pdf_raw(const Matrix& u) override;

  Eigen::VectorXd pdf(const Matrix& u) override;

  Eigen::VectorXd cdf(const Matrix& u) override;

  Eigen::VectorXd hfunc1_raw(const Matrix& u) override;

  Eigen::VectorXd hfunc2_raw(const Matrix& u) override;

  Eigen::VectorXd hfunc1(const Matrix& u) override;

  Eigen::VectorXd hfunc2(const Matrix& u) override;

  Eigen::VectorXd hinv1_raw(const Matrix& u) override;

  Eigen::VectorXd hinv2_raw(const Matrix& u) override;

  double get_npars() const override;

  void set_npars(const double& npars) override;

  Matrix get_parameters() const override;

  Matrix get_parameters_lower_bounds() const override;

  Matrix get_parameters_upper_bounds() const override;

  void set_parameters(const Matrix& parameters) override;

  double parameters_to_tau(const Matrix& parameters) override;

  void flip() override;

  Matrix tau_to_parameters(const double& tau) override;

  Eigen::VectorXd make_normal_grid(size_t m = 30);

  std::shared_ptr<tools_interpolation::InterpolationGrid> interp_grid_;
  double npars_;
};
}

#include <vinecopulib/bicop/implementation/kernel.ipp>
