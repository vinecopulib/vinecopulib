// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

namespace vinecopulib
{
    class IndepBicop : public ParBicop
    {
    public:
        // constructor
        IndepBicop();
        IndepBicop(const Eigen::VectorXd& parameters);
        IndepBicop(const Eigen::VectorXd& parameters, const int& rotation);

        // PDF
        Eigen::VectorXd pdf_default(const Eigen::MatrixXd& u);

        // hfunctions and their inverses
        Eigen::VectorXd hfunc1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u);

        Eigen::VectorXd tau_to_parameters(const double &);
        double parameters_to_tau(const Eigen::VectorXd &);

        void flip();

    private:
        Eigen::VectorXd get_start_parameters(const double tau);
    };
}
