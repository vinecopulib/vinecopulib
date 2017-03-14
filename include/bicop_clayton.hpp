// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_archimedean.hpp"

class ClaytonBicop : public ArchimedeanBicop {

public:
    // constructor
    ClaytonBicop();
    ClaytonBicop(const VecXd& parameters);
    ClaytonBicop(const VecXd& parameters, const int& rotation);

    // generator, its inverse and derivatives for the archimedean copula
    VecXd generator(const VecXd& u);
    VecXd generator_inv(const VecXd& u);
    VecXd generator_derivative(const VecXd& u);
    VecXd generator_derivative2(const VecXd& u);

    // inverse hfunction
    VecXd hinv1_default(const MatXd& u);

    // link between Kendall's tau and the par_bicop parameter
    VecXd tau_to_parameters(const double& tau);
    double parameters_to_tau(const VecXd& parameters);

private:
    VecXd get_start_parameters(const double tau);
};
