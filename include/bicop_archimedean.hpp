// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

namespace vinecopulib
{
    class ArchimedeanBicop : public ParBicop {

    public:
        // pdf, hfunctions and inverses
        VectorXd pdf_default(const MatrixXd& u);
        VectorXd hfunc1_default(const MatrixXd& u);
        VectorXd hfunc2_default(const MatrixXd& u);
        VectorXd hinv1_default(const MatrixXd& u);
        VectorXd hinv2_default(const MatrixXd& u);

        // generator, its inverse and derivatives
        virtual VectorXd generator(const VectorXd& u) = 0;
        virtual VectorXd generator_inv(const VectorXd& u) = 0;
        virtual VectorXd generator_derivative(const VectorXd& u) = 0;
        virtual VectorXd generator_derivative2(const VectorXd& u) = 0;

        void flip();

    private:
        VectorXd get_start_parameters(const double tau);
    };
}

