// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/circular.hpp>

namespace vinecopulib {

//! @brief The circular-linear copula with cubic sections in the linear
//! variable.
//!
//! The circular variable is the first argument internally; when `var_types`
//! puts it second, the leaves are evaluated on swapped columns with the two
//! h-functions (and their inverses) exchanged.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class CubicSectionsBicop : public CircularBicop
{
public:
  CubicSectionsBicop();

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

private:
  using Leaf = double (CubicSectionsBicop::*)(
    double,
    double,
    const Eigen::Ref<const Eigen::VectorXd>&) const;

  // evaluates `leaf` row by row with the circular variable first; when it is
  // the second argument, the columns are swapped and `swapped` is used instead
  Eigen::VectorXd eval(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters,
                       Leaf leaf,
                       Leaf swapped) const;

  // the leaves with the circular variable `u` as the first argument
  double pdf_at(double u,
                double v,
                const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double cdf_at(double u,
                double v,
                const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hfunc1_at(double u,
                   double v,
                   const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hfunc2_at(double u,
                   double v,
                   const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hinv1_at(double u,
                  double w,
                  const Eigen::Ref<const Eigen::VectorXd>& par) const;
  double hinv2_at(double w,
                  double v,
                  const Eigen::Ref<const Eigen::VectorXd>& par) const;

  // the section polynomial p(v) and its antiderivative P(v), P(0) = P(1) = 0
  static double p_of(double a, double b, double v);
  static double P_of(double a, double b, double v);
};
}

#include <vinecopulib/bicop/implementation/cubic_sections.ipp>
