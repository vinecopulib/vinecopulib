// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_archimedean.hpp"

namespace vinecopulib
{
    class JoeBicop : public ArchimedeanBicop
    {
    public:
        // constructor
        JoeBicop();
        JoeBicop(const Eigen::VectorXd& parameters);
        JoeBicop(const Eigen::VectorXd& parameters, const int& rotation);

        // generator, its inverse and derivatives for the archimedean copula
        double generator(const double& u);
        double generator_inv(const double& u);
        double generator_derivative(const double& u);
        double generator_derivative2(const double& u);

        // inverse hfunction
        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u);

        // link between Kendall's tau and the par_bicop parameter
        Eigen::VectorXd tau_to_parameters(const double& tau);
        double parameters_to_tau(const Eigen::VectorXd& par);

    private:
        Eigen::VectorXd get_start_parameters(const double tau);
    };
}

double qcondjoe(double* q, double* u, double* de);
