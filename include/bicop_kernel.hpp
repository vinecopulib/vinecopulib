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
    class KernelBicop : public Bicop {
    public:
        KernelBicop();

        VectorXd pdf_default(const MatrixXd& u);
        VectorXd hfunc1_default(const MatrixXd& u);
        VectorXd hfunc2_default(const MatrixXd& u);
        VectorXd hinv1_default(const MatrixXd& u);
        VectorXd hinv2_default(const MatrixXd& u);

        double parameters_to_tau(const VectorXd &);
        double calculate_npars();

        void flip();

    protected:
        InterpolationGrid interp_grid_;
        double npars_;
    };
}
