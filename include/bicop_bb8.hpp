// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_archimedean.hpp"

namespace vinecopulib
{
    class Bb8Bicop : public ArchimedeanBicop {

    public:
        // constructor
        Bb8Bicop();
        Bb8Bicop(const VecXd& parameters);
        Bb8Bicop(const VecXd& parameters, const int& rotation);

        // generator, its inverse and derivatives for the archimedean copula
        VecXd generator(const VecXd& u);
        VecXd generator_inv(const VecXd& u);
        VecXd generator_derivative(const VecXd& u);
        VecXd generator_derivative2(const VecXd& u);

        // link between Kendall's tau and the par_bicop parameter
        double parameters_to_tau(const VecXd& par);
    };
}
