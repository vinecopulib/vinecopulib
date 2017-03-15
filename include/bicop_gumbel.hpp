// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_archimedean.hpp"

namespace vinecopulib
{
    class GumbelBicop : public ArchimedeanBicop {
    public:
        // constructor
        GumbelBicop();
        GumbelBicop(const VectorXd& parameters);
        GumbelBicop(const VectorXd& parameters, const int& rotation);

        // generator, its inverse and derivatives for the archimedean copula
        VectorXd generator(const VectorXd& u);
        VectorXd generator_inv(const VectorXd& u);
        VectorXd generator_derivative(const VectorXd& u);
        VectorXd generator_derivative2(const VectorXd& u);

        // inverse hfunction
        VectorXd hinv1_default(const MatrixXd& u);

        // link between Kendall's tau and the par_bicop parameter
        VectorXd tau_to_parameters(const double& tau);
        double parameters_to_tau(const VectorXd& parameters);

    private:
        VectorXd get_start_parameters(const double tau);
    };
}

double qcondgum(double* q, double* u, double* de);
