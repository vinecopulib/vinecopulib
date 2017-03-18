// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_archimedean.hpp"

namespace vinecopulib
{
    class ClaytonBicop : public ArchimedeanBicop
    {
    public:
        ClaytonBicop();

        Eigen::VectorXd generator(const Eigen::VectorXd& u);
        Eigen::VectorXd generator_inv(const Eigen::VectorXd& u);
        Eigen::VectorXd generator_derivative(const Eigen::VectorXd& u);
        Eigen::VectorXd generator_derivative2(const Eigen::VectorXd& u);

        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u);

        std::vector<Eigen::MatrixXd> tau_to_parameters(const double& tau);
        double parameters_to_tau(const std::vector<Eigen::MatrixXd>& parameters);

    private:
        std::vector<Eigen::MatrixXd> get_start_parameters(const double tau);
    };
}
