// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "misc/interpolation_grid.hpp"
#include "bicop/class.hpp"

namespace vinecopulib
{
    class KernelBicop : public Bicop
    {
    public:
        KernelBicop();

    private:
        Eigen::VectorXd pdf_default(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hfunc1_default(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hfunc2_default(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hinv1_default(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hinv2_default(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        double parameters_to_tau(const Eigen::VectorXd &);
        Eigen::MatrixXd tau_to_parameters_default(const double& tau);
        double calculate_npars();

        void flip();

    protected:
        InterpolationGrid interp_grid_;
        double npars_;
    };
}
