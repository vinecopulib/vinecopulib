// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop/parametric.hpp"

namespace vinecopulib
{
    class ArchimedeanBicop : public ParBicop
    {
    private:
        // pdf, hfunctions and inverses
        Eigen::VectorXd pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);

        // generator, its inverse and derivatives
        virtual double generator(const double& u) = 0;
        virtual double generator_inv(const double& u) = 0;
        virtual double generator_derivative(const double& u) = 0;
        virtual double generator_derivative2(const double& u) = 0;

        Eigen::VectorXd get_start_parameters(const double tau);
    protected:
        double flip_tau(double tau);
    };
}
