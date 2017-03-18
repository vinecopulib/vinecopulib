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

        // PDF
        Eigen::VectorXd pdf_default(const Eigen::MatrixXd& u);

        // hfunctions and their inverses
        Eigen::VectorXd hfunc1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u);

        std::vector<Eigen::MatrixXd> tau_to_parameters(const double &);
        double parameters_to_tau(const std::vector<Eigen::MatrixXd>&);

        void flip();

    private:
        std::vector<Eigen::MatrixXd> get_start_parameters(const double tau);
    };
}
