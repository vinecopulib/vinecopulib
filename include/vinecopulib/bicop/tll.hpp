// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/kernel.hpp>

namespace vinecopulib {
//! @brief The transformation local likelihood estimator (`tll`).
class TllBicop : public KernelBicop
{
public:
  TllBicop();

private:
  static Eigen::VectorXd gaussian_kernel_2d(const Eigen::MatrixXd& x);

  Eigen::Matrix2d select_bandwidth(const Eigen::MatrixXd& x,
                                   const std::string& method,
                                   const Eigen::VectorXd& weights);

  Eigen::MatrixXd fit_local_likelihood(const Eigen::MatrixXd& x,
                                       const Eigen::MatrixXd& x_data,
                                       const Eigen::Matrix2d& B,
                                       const std::string& method,
                                       const Eigen::VectorXd& weights);

  double calculate_infl(const size_t& n,
                        const double& f0,
                        const Eigen::Vector2d& b,
                        const Eigen::Matrix2d& B,
                        const double& det_irB,
                        const Eigen::Matrix2d& S,
                        const std::string& method,
                        const double& weight);

  void fit(const Eigen::MatrixXd& data,
           std::string method,
           double mult,
           size_t grid_size,
           const Eigen::VectorXd& weights) override;

  // the estimator for a pair with a circular variable: a product kernel that
  // is von Mises on a circular axis and Gaussian on the probit scale of a
  // linear one, with a local log-linear correction per axis
  //! @brief One axis of the mixed-geometry estimator.
  struct Axis
  {
    bool circular;
    double scale;             // von Mises concentration or Gaussian sd
    Eigen::VectorXd x;        // transformed observations
    Eigen::VectorXd grid;     // transformed knots
    Eigen::VectorXd knots;    // knots on the copula scale
    Eigen::VectorXd jacobian; // density of the transform at the knots
  };

  void fit_mixed(const Eigen::MatrixXd& data,
                 const std::string& method,
                 double mult,
                 size_t grid_size,
                 const Eigen::VectorXd& weights);

  static Axis make_axis(const std::string& var_type,
                        const Eigen::VectorXd& u,
                        size_t grid_size);

  static double select_bandwidth_mixed(const Axis& axis,
                                       size_t n,
                                       const std::string& method,
                                       double dependence);

  static std::pair<double, double> local_fit_mixed(
    const std::vector<Axis>& axes,
    const std::array<Eigen::Index, 2>& knot,
    const std::string& method,
    const Eigen::VectorXd& weights,
    Eigen::VectorXd& kernels);
};
}

#include <vinecopulib/bicop/implementation/tll.ipp>
