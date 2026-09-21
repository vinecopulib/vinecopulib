// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/binding.hpp>

namespace vinecopulib {
//! @brief The wrapped Cauchy circula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! The binding density is the wrapped Cauchy
//! \f$ g(\theta) = (1 - \rho^2) / (2\pi(1 + \rho^2 - 2\rho\cos\theta)) \f$
//! with concentration \f$ \rho \in [0, 0.99] \f$; its lifted CDF and inverse
//! are closed form.
//!
//! @literature
//! Jones, M. C., Pewsey, A., and Kato, S. (2015). On a class of circulas:
//! copulas for circular distributions. Annals of the Institute of Statistical
//! Mathematics, 67, 843–862.
class WrappedCauchyBicop : public BindingBicop
{
public:
  WrappedCauchyBicop();

private:
  double g(double theta, double concentration) const override;

  double lifted_cdf(double theta, double concentration) const override;

  double lifted_cdf_inverse(double w, double concentration) const override;

  double concentration_from_resultant(double rbar) const override;
};
}

#include <vinecopulib/bicop/implementation/wrapped_cauchy.ipp>
