// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vector>
#include <vinecopulib/bicop/circular.hpp>

namespace vinecopulib {

//! @brief An abstract class for the binding-density circulas.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! The copula density is \f$ c(u, v) = 2\pi g(2\pi(v - u) - \mu) \f$ for a
//! circular density \f$ g \f$ symmetric about zero with concentration
//! parameter `parameters(0)` and phase \f$ \mu = \f$ `parameters(1)`. The
//! opposite orientation is the 90-degree rotation. Every leaf is expressed
//! through the lifted distribution function \f$ \tilde G(\theta) =
//! \int_0^\theta g \f$ and its inverse on the whole real line.
//!
//! @literature
//! Jones, M. C., Pewsey, A., and Kato, S. (2015). On a class of circulas:
//! copulas for circular distributions. Annals of the Institute of Statistical
//! Mathematics, 67, 843–862.
class BindingBicop : public CircularBicop
{
protected:
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  void flip() override;

  Eigen::VectorXd moment_start(const Eigen::MatrixXd& data,
                               const Eigen::VectorXd& weights) override;

  // the circular density, its lifted CDF, and the inverse of the lifted CDF,
  // all as functions of the angle and the concentration
  virtual double g(double theta, double concentration) const = 0;

  virtual double lifted_cdf(double theta, double concentration) const = 0;

  virtual double lifted_cdf_inverse(double w, double concentration) const = 0;

  // the concentration whose mean resultant length is `rbar`
  virtual double concentration_from_resultant(double rbar) const = 0;

  // the cosine coefficients rho_j, j >= 1, of the circular density
  // g(theta) = (1 + 2 sum_j rho_j cos(j theta)) / (2 pi), truncated where
  // they no longer matter; they give the CDF of the copula in closed form
  virtual std::vector<double> fourier_coefficients(
    double concentration) const = 0;
};
}

#include <vinecopulib/bicop/implementation/binding.ipp>
