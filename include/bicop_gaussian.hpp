// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_elliptical.hpp"

namespace vinecopulib
{
    class GaussianBicop : public EllipticalBicop {

    public:
        // constructor
        GaussianBicop();
        GaussianBicop(const VecXd& parameters);
        GaussianBicop(const VecXd& parameters, const int& rotation);

        // PDF
        VecXd pdf_default(const MatXd& u);

        // hfunction
        VecXd hfunc1_default(const MatXd& u);

        // inverse hfunction
        VecXd hinv1_default(const MatXd& u);

    private:
        VecXd get_start_parameters(const double tau);
    };
}
