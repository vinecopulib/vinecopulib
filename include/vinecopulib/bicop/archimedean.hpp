// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <cmath>
#include <limits>
#include <vinecopulib/bicop/parametric.hpp>

namespace vinecopulib {
//! @brief An abstract class for Archimedean copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class ArchimedeanBicop : public ParBicop
{
protected:
  // generic generator-based cdf; families override with a vectorized
  // closed form and defer to this for the per-row-parameters case
  Eigen::VectorXd cdf(const tools_eigen::ConstMatRef& u,
                      const tools_eigen::ConstMatRef& parameters);

  //! Shared scaffolding for the closed-form h-functions: the
  //! broadcast/per-row branch plus the common NaN + clamp postprocess (NaN
  //! input -> NaN, a numerically failed evaluation -> `uother`, else
  //! `min(h, 1)`). Each family supplies only its math kernels:
  //!   - `broadcast(uc, uo)`: array kernel over the whole column (single
  //!     broadcast parameter set), returning an `Eigen::ArrayXd`;
  //!   - `scalar(i, uc, uo)`: per-row kernel reading `parameters.row(i)`,
  //!     returning a `double`.
  //! Zero-cost vs. hand-written per family: the broadcast kernel still yields
  //! one fused expression and the scalar kernel inlines through the template.
  template<typename Broadcast, typename Scalar>
  static Eigen::VectorXd apply_closed_form_h(
    const Eigen::Ref<const Eigen::VectorXd>& ucond,
    const Eigen::Ref<const Eigen::VectorXd>& uother,
    const tools_eigen::ConstMatRef& parameters,
    Broadcast&& broadcast,
    Scalar&& scalar)
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (parameters.rows() == 1) {
      const auto uc = ucond.array();
      const auto uo = uother.array();
      Eigen::ArrayXd h = broadcast(uc, uo);
      return ((uc.isNaN() || uo.isNaN())
                .select(nan, h.isNaN().select(uo, h.min(1.0))))
        .matrix();
    }
    const Eigen::Index n = ucond.size();
    Eigen::VectorXd out(n);
    for (Eigen::Index i = 0; i < n; ++i) {
      const double uc = ucond(i);
      const double uo = uother(i);
      if ((std::isnan)(uc) || (std::isnan)(uo)) {
        out(i) = nan;
        continue;
      }
      const double h = scalar(i, uc, uo);
      out(i) = (std::isnan)(h) ? uo : std::min(h, 1.0);
    }
    return out;
  }

private:
  // hfunctions and inverses (`parameters` is m x p, m in {1, n}; a single
  // row is broadcast to all observations)
  Eigen::VectorXd hfunc1_raw(const tools_eigen::ConstMatRef& u,
                             const tools_eigen::ConstMatRef& parameters);

  Eigen::VectorXd hfunc2_raw(const tools_eigen::ConstMatRef& u,
                             const tools_eigen::ConstMatRef& parameters);

  Eigen::VectorXd hinv1_raw(const tools_eigen::ConstMatRef& u,
                            const tools_eigen::ConstMatRef& parameters);

  Eigen::VectorXd hinv2_raw(const tools_eigen::ConstMatRef& u,
                            const tools_eigen::ConstMatRef& parameters);

  // Archimedean copulas are exchangeable: the second h-function derivatives
  // are the first ones at swapped arguments/selectors
  Eigen::VectorXd hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv);

  Eigen::VectorXd hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv);

  // generator, its inverse and derivative; `parameters` is a single parameter
  // set (a p x 1 column)
  virtual double generator(
    const double& u,
    const Eigen::Ref<const Eigen::VectorXd>& parameters) = 0;

  virtual double generator_inv(
    const double& u,
    const Eigen::Ref<const Eigen::VectorXd>& parameters) = 0;

  virtual double generator_derivative(
    const double& u,
    const Eigen::Ref<const Eigen::VectorXd>& parameters) = 0;

  // virtual double generator_derivative2(const double &u) = 0;

  Eigen::VectorXd get_start_parameters(const double tau);
};
}

#include <vinecopulib/bicop/implementation/archimedean.ipp>
