// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

class ArchimedeanBicop : public ParBicop {

public:
    // pdf, hfunctions and inverses
    VecXd pdf_default(const MatXd& u);
    VecXd hfunc1_default(const MatXd& u);
    VecXd hfunc2_default(const MatXd& u);
    VecXd hinv1_default(const MatXd& u);
    VecXd hinv2_default(const MatXd& u);

    // generator, its inverse and derivatives
    virtual VecXd generator(const VecXd& u) = 0;
    virtual VecXd generator_inv(const VecXd& u) = 0;
    virtual VecXd generator_derivative(const VecXd& u) = 0;
    virtual VecXd generator_derivative2(const VecXd& u) = 0;

    void flip();

private:
    VecXd get_start_parameters(const double tau);
};
