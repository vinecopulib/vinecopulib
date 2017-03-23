// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop/archimedean.hpp"

namespace vinecopulib
{
    class ClaytonBicop : public ArchimedeanBicop
    {
    public:
        // constructor
        ClaytonBicop();

    private:
        // generator, its inverse and derivatives for the archimedean copula
        double generator(const double& u);
        double generator_inv(const double& u);
        double generator_derivative(const double& u);
        double generator_derivative2(const double& u);

        // inverse hfunction
        Eigen::VectorXd hinv1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        // link between Kendall's tau and the par_bicop parameter
        Eigen::MatrixXd tau_to_parameters(const double& tau);
        double parameters_to_tau(const Eigen::VectorXd& parameters);

        Eigen::VectorXd get_start_parameters(const double tau);
    };
}
