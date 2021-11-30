// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/kernel.hpp>

namespace vinecopulib {
//! @brief The transformation local-constant likelihood estimator.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Nagler, Thomas. *kdecopula: An R Package for the Kernel Estimation of
//! Copula Densities*. arXiv:1603.04229 [stat.CO], 2016
class TllBicop : public KernelBicop
{
public:
  TllBicop();

private:
  Vector gaussian_kernel_2d(const Matrix& x);

  Eigen::Matrix2d select_bandwidth(const Matrix& x,
                                   std::string method,
                                   const Vector& weights);

  Matrix fit_local_likelihood(const Matrix& x,
                              const Matrix& x_data,
                              const Eigen::Matrix2d& B,
                              std::string method,
                              const Vector& weights);

  double calculate_infl(const size_t& n,
                        const double& f0,
                        const Eigen::Vector2d& b,
                        const Eigen::Matrix2d& B,
                        const double& det_irB,
                        const Eigen::Matrix2d& S,
                        const std::string& method,
                        const double& weight);

  void fit(const Matrix& data,
           std::string method,
           double mult,
           const Vector& weights);
};
}

#include <vinecopulib/bicop/implementation/tll.ipp>
