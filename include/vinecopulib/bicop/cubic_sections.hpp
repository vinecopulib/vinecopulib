// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/sections.hpp>

namespace vinecopulib {
//! @brief The circular-linear copula with cubic sections.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! \f$ c(u, v) = 1 + \cos(2\pi u - \mu)\bigl[a (1 - v)(1 - 3v) + b\, v (2 -
//! 3v)\bigr] \f$ with amplitudes \f$ a \in [0, 1] \f$ and \f$ b \in [-1, 1]
//! \f$ and phase \f$ \mu \f$. The amplitude \f$ a \f$ controls the direction
//! preference near \f$ v = 0 \f$ and \f$ b \f$ the one near \f$ v = 1 \f$.
//!
//! @literature
//! Hodel, F. H. and Fieberg, J. R. (2022). Circular-linear copulae for animal
//! movement data. Methods in Ecology and Evolution, 13, 1001–1013.
class CubicSectionsBicop : public SectionsBicop
{
public:
  CubicSectionsBicop();

private:
  std::pair<double, double> amplitudes(
    const Eigen::Ref<const Eigen::VectorXd>& parameters) const override;

  Eigen::VectorXd parameters_from_moments(double a,
                                          double b,
                                          double mu) const override;
};
}

#include <vinecopulib/bicop/implementation/cubic_sections.ipp>
