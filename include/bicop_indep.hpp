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
        IndepBicop(const VecXd& parameters);
        IndepBicop(const VecXd& parameters, const int& rotation);

        // PDF
        VecXd pdf_default(const MatXd& u);

        // hfunctions: the conditioning variable is put second
        VecXd hfunc1_default(const MatXd& u);
        VecXd hfunc2_default(const MatXd& u);
        VecXd hinv1_default(const MatXd& u);
        VecXd hinv2_default(const MatXd& u);

        VecXd tau_to_parameters(const double &);
        double parameters_to_tau(const VecXd &);

        void flip();

    private:
        VecXd get_start_parameters(const double tau);
    };
}
