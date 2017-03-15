// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

namespace vinecopulib
{
    class EllipticalBicop : public ParBicop  {
    public:
        // hfunction and its inverse
        VecXd hfunc2_default(const MatXd& u);
        VecXd hinv2_default(const MatXd& u);

        // link between Kendall's tau and the par_bicop parameter
        VecXd tau_to_parameters(const double& tau);
        double parameters_to_tau(const VecXd& parameters);

        void flip();
    };
}
