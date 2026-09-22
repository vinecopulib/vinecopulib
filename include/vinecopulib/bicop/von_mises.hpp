// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/binding.hpp>

namespace vinecopulib {
//! @brief The von Mises circula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! The binding density is the von Mises
//! \f$ g(\theta) = e^{\kappa\cos\theta} / (2\pi I_0(\kappa)) \f$ with
//! concentration \f$ \kappa \in [0, 100] \f$. Its lifted CDF is a Fourier
//! series in Bessel-function ratios and the inverse a bracketed root solve.
//!
//! @literature
//! Jones, M. C., Pewsey, A., and Kato, S. (2015). On a class of circulas:
//! copulas for circular distributions. Annals of the Institute of Statistical
//! Mathematics, 67, 843–862.
class VonMisesBicop : public BindingBicop
{
public:
  VonMisesBicop();

private:
  double g(double theta, double concentration) const override;

  double lifted_cdf(double theta, double concentration) const override;

  double lifted_cdf_inverse(double w, double concentration) const override;

  double concentration_from_resultant(double rbar) const override;

  std::vector<double> fourier_coefficients(double concentration) const override;

  //! @brief The Bessel-function ratios of the von Mises series for one
  //! concentration, kept between consecutive calls with the same value.
  struct Series
  {
    //! the concentration the ratios belong to; negative when unset
    double kappa = -1.0;
    //! the ratios \f$ I_k(\kappa) / I_0(\kappa) \f$ for \f$ k \ge 1 \f$
    std::vector<double> ratios;
  };
  static const std::vector<double>& ratios_for(Series& series, double kappa);
};
}

#include <vinecopulib/bicop/implementation/von_mises.ipp>
