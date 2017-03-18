// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

namespace vinecopulib
{
    class ArchimedeanBicop : public ParBicop
    {
    public:
        Eigen::VectorXd pdf_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u);

        virtual Eigen::VectorXd generator(const Eigen::VectorXd& u) = 0;
        virtual Eigen::VectorXd generator_inv(const Eigen::VectorXd& u) = 0;
        virtual Eigen::VectorXd generator_derivative(const Eigen::VectorXd& u) = 0;
        virtual Eigen::VectorXd generator_derivative2(const Eigen::VectorXd& u) = 0;

        void flip();

    protected:
        double flip_tau(double tau);
    private:
        std::vector<Eigen::MatrixXd> get_start_parameters(const double tau);
    };
}
