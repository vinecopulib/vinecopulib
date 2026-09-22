// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
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

class KernelBicop : public AbstractBicop
{
public:
  KernelBicop();

  void set_var_types(const std::vector<std::string>& var_types) override;

  std::vector<Eigen::VectorXd> get_grid_knots() const override;

  void set_grid(const std::vector<Eigen::VectorXd>& knots,
                const Eigen::MatrixXd& values) override;

protected:
  // evaluation leaves; kernel estimators store an interpolation grid rather
  // than a per-row parameter vector, so they ignore `parameters`
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters) override;

  // the interpolant is piecewise bilinear, so these two are exact sums of
  // nonnegative weights against a nonnegative grid rather than differences
  double rect_prob(double a1,
                   double b1,
                   double a2,
                   double b2,
                   const Eigen::MatrixXd& parameters) override;

  double cond_interval_prob(double u_cond,
                            double lo,
                            double hi,
                            size_t cond_var,
                            const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  double get_npars() const override;

  void set_npars(const double& npars) override;

  Eigen::MatrixXd get_parameters() const override;

  Eigen::MatrixXd get_parameters_lower_bounds() const override;

  Eigen::MatrixXd get_parameters_upper_bounds() const override;

  void set_parameters(const Eigen::MatrixXd& parameters) override;

  double parameters_to_tau(const Eigen::MatrixXd& parameters) override;

  void flip() override;

  Eigen::MatrixXd tau_to_parameters(const double& tau) override;

  // the default knots of an axis: equally spaced on the normal scale for a
  // linear axis, equally spaced on the unit interval for a circular one
  static Eigen::VectorXd make_grid_points(const std::string& var_type,
                                          size_t m);

  Eigen::VectorXd make_normal_grid(size_t m = 30);

  std::shared_ptr<tools_interpolation::InterpolationGrid> interp_grid_;
  double npars_;
};
}

#include <vinecopulib/bicop/implementation/kernel.ipp>
