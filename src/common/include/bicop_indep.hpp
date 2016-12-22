/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecopulib.

    vinecoplib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecoplib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VINECOPLIB_BICOP_INDEP_HPP
#define VINECOPLIB_BICOP_INDEP_HPP

#include "bicop_parametric.hpp"

class IndepBicop : public ParBicop {

public:
    // constructor
    IndepBicop();

    // PDF
    VecXd pdf(const MatXd& u);

    // hfunctions: the conditioning variable is put second
    VecXd hfunc1(const MatXd& u);
    VecXd hfunc2(const MatXd& u);

    VecXd tau_to_par(const double __attribute__((unused))& tau);
    double par_to_tau(const VecXd __attribute__((unused))& parameters);

    VecXd hinv1(const MatXd& u);
    VecXd hinv2(const MatXd& u);
};

#endif
