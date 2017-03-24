// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "misc/interpolation_grid.hpp"
#include "bicop/abstract.hpp"

namespace vinecopulib
{
    class KernelBicop : public AbstractBicop
    {
    public:
        KernelBicop();

    protected:
        Eigen::VectorXd pdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hfunc1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hfunc2(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hinv1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hinv2(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        double parameters_to_tau(const Eigen::VectorXd &);
        Eigen::MatrixXd tau_to_parameters(const double& tau);
        double calculate_npars();

        void flip();

        InterpolationGrid interp_grid_;
        double npars_;
    };
}
