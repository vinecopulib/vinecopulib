// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop/elliptical.hpp"

namespace vinecopulib
{
    class StudentBicop : public EllipticalBicop
    {
    public:
        // constructor
        StudentBicop();

    private:
        // PDF
        Eigen::VectorXd pdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        // hfunction
        Eigen::VectorXd hfunc1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        // inverse hfunction
        Eigen::VectorXd hinv1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        Eigen::MatrixXd tau_to_parameters(const double& tau);

        Eigen::VectorXd get_start_parameters(const double tau);
    };
}
