// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/circular.hpp>

namespace vinecopulib {

//! @brief An abstract class for the circular-linear copulas with polynomial
//! sections in the linear variable.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! With the circular variable \f$ u \f$ and the linear variable \f$ v \f$, the
//! density is \f$ c(u, v) = 1 + \cos(2\pi u - \mu)\, p(v) \f$ with
//! \f$ p(v) = a (1 - v)(1 - 3v) + b\, v (2 - 3v) \f$; the quadratic family is
//! the case \f$ a = b \f$. The leaves read which argument is circular from
//! the variable types, so both argument orders are supported with the same
//! parameters.
//!
//! @literature
//! Hodel, F. H. and Fieberg, J. R. (2022). Circular-linear copulae for animal
//! movement data. Methods in Ecology and Evolution, 13, 1001–1013.
class SectionsBicop : public CircularBicop
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

  double parameters_to_tau(const Eigen::MatrixXd& parameters) override;

  void flip() override;

  Eigen::VectorXd moment_start(const Eigen::MatrixXd& data,
                               const Eigen::VectorXd& weights) override;

  // the amplitudes (a, b) of the two ends of the linear variable
  virtual std::pair<double, double> amplitudes(
    const Eigen::Ref<const Eigen::VectorXd>& parameters) const = 0;

  // the parameters implied by the moment estimates of (a, b) and the phase
  virtual Eigen::VectorXd parameters_from_moments(double a,
                                                  double b,
                                                  double mu) const = 0;

private:
  // whether the circular variable is the second argument
  bool circular_second() const;

  // the leaves with the circular variable as the first argument
  double pdf_internal(double u,
                      double v,
                      const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double cdf_internal(double u,
                      double v,
                      const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hfunc1_internal(double u,
                         double v,
                         const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hfunc2_internal(double u,
                         double v,
                         const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hinv1_internal(double u,
                        double w,
                        const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hinv2_internal(double w,
                        double v,
                        const Eigen::Ref<const Eigen::VectorXd>& par) const;

  static double p_of(double a, double b, double v);
  static double P_of(double a, double b, double v);
};
}

#include <vinecopulib/bicop/implementation/sections.ipp>
