// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "interpolation_grid.hpp"
#include "bicop_class.hpp"

class KernelBicop : public Bicop {
public:
    KernelBicop();

    VecXd pdf_default(const MatXd& u);
    VecXd hfunc1_default(const MatXd& u);
    VecXd hfunc2_default(const MatXd& u);
    VecXd hinv1_default(const MatXd& u);
    VecXd hinv2_default(const MatXd& u);

    double parameters_to_tau(const VecXd &);
    double calculate_npars();

    void flip();

protected:
    InterpolationGrid interp_grid_;
    double npars_;
};
