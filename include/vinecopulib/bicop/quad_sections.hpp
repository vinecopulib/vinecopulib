// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/sections.hpp>

namespace vinecopulib {
//! @brief The circular-linear copula with quadratic sections.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! \f$ c(u, v) = 1 + a \cos(2\pi u - \mu)(1 - 2v) \f$ with amplitude
//! \f$ a \in [0, 1] \f$ and phase \f$ \mu \f$; the circular-linear analog of
//! the Farlie-Gumbel-Morgenstern copula.
//!
//! @literature
//! Hodel, F. H. and Fieberg, J. R. (2022). Circular-linear copulae for animal
//! movement data. Methods in Ecology and Evolution, 13, 1001–1013.
class QuadSectionsBicop : public SectionsBicop
{
public:
  QuadSectionsBicop();

private:
  std::pair<double, double> amplitudes(
    const Eigen::Ref<const Eigen::VectorXd>& parameters) const override;

  Eigen::VectorXd parameters_from_moments(double a,
                                          double b,
                                          double mu) const override;
};
}

#include <vinecopulib/bicop/implementation/quad_sections.ipp>
