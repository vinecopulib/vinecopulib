// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

namespace vinecopulib
{
    class IndepBicop : public ParBicop {
    public:
        // constructor
        IndepBicop();
        IndepBicop(const VectorXd& parameters);
        IndepBicop(const VectorXd& parameters, const int& rotation);

        // PDF
        VectorXd pdf_default(const MatrixXd& u);

        // hfunctions: the conditioning variable is put second
        VectorXd hfunc1_default(const MatrixXd& u);
        VectorXd hfunc2_default(const MatrixXd& u);
        VectorXd hinv1_default(const MatrixXd& u);
        VectorXd hinv2_default(const MatrixXd& u);

        VectorXd tau_to_parameters(const double &);
        double parameters_to_tau(const VectorXd &);

        void flip();

    private:
        VectorXd get_start_parameters(const double tau);
    };
}
