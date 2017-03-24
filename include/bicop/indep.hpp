// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop/parametric.hpp"

namespace vinecopulib
{
    class IndepBicop : public ParBicop
    {
    public:
        // constructor
        IndepBicop();

    private:
        // PDF
        Eigen::VectorXd pdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        // hfunctions and their inverses
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

        Eigen::MatrixXd tau_to_parameters(const double &);
        double parameters_to_tau(const Eigen::VectorXd &);

        void flip();

        Eigen::VectorXd get_start_parameters(const double tau);
    };
}
