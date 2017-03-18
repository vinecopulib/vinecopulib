// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "interpolation_grid.hpp"
#include "bicop_class.hpp"

namespace vinecopulib
{
    class KernelBicop : public Bicop
    {
    public:
        KernelBicop();

        Eigen::VectorXd pdf_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u);

        double parameters_to_tau(const std::vector<Eigen::MatrixXd>&);
        double calculate_npars();

        void flip();

    protected:
        InterpolationGrid interp_grid_;
        double npars_;
    };
}
