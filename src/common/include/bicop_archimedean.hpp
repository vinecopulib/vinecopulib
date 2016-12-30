/*
    Copyright 2016 Thibault Vatter, Thomas Nagler

    This file is part of vinecopulib.

    vinecopulib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecopulib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VINECOPULIB_BICOP_ARCHIMEDEAN_HPP
#define VINECOPULIB_BICOP_ARCHIMEDEAN_HPP

#include "bicop_parametric.hpp"

class ArchimedeanBicop : public ParBicop {

public:
    // PDF
    VecXd pdf_default(const MatXd& u);

    // hfunctions: the conditioning variable is put second
    VecXd hfunc1_default(const MatXd& u);
    VecXd hfunc2_default(const MatXd& u);
    VecXd hinv1_default(const MatXd& u);
    VecXd hinv2_default(const MatXd& u);

    // generator, its inverse and derivatives for the archimedean copula
    virtual VecXd generator(const VecXd& u) = 0;
    virtual VecXd generator_inv(const VecXd& u) = 0;
    virtual VecXd generator_derivative(const VecXd& u) = 0;
    virtual VecXd generator_derivative2(const VecXd& u) = 0;

    virtual VecXd hinv(const MatXd& u) = 0;
};


#endif
