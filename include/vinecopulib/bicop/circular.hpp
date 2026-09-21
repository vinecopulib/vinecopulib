// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/misc/tools_circular.hpp>

namespace vinecopulib {

//! @brief An abstract class for the parametric families with a circular
//! variable and a phase parameter.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! The last parameter of every derived family is a phase in radians. It is
//! unbounded for the optimizer and reduced to \f$ [-\pi, \pi) \f$ after a
//! fit. Fitting starts from moment estimates instead of Kendall's tau, which
//! carries no information about circular association, and Kendall's tau of
//! the fitted model is computed by numerical integration.
class CircularBicop : public ParBicop
{
protected:
  void fit(const Eigen::MatrixXd& data,
           std::string method,
           double,
           size_t,
           const Eigen::VectorXd& weights) override;

  double parameters_to_tau(const Eigen::MatrixXd& parameters) override;

  Eigen::MatrixXd parameters_to_taildep(
    const Eigen::MatrixXd& parameters) override;

  Eigen::MatrixXd tau_to_parameters(const double& tau) override;

  Eigen::VectorXd get_start_parameters(const double tau) override;

  // moment estimates of all parameters from (rotated, trimmed) copula data;
  // the starting values of the maximum-likelihood fit
  virtual Eigen::VectorXd moment_start(const Eigen::MatrixXd& data,
                                       const Eigen::VectorXd& weights) = 0;

  // the phase is the last parameter
  Eigen::Index phase_index() const;
};
}

#include <vinecopulib/bicop/implementation/circular.ipp>
