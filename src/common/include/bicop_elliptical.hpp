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

#ifndef VINECOPULIB_BICOP_ELLIPTICAL_HPP
#define VINECOPULIB_BICOP_ELLIPTICAL_HPP

#include "bicop_parametric.hpp"

class EllipticalBicop : public ParBicop  {
public:
    // copula density
    virtual VecXd pdf_default(const MatXd& u) = 0;

    // hfunction and its inverse
    virtual VecXd hfunc1_default(const MatXd& u) = 0;
    VecXd hfunc2_default(const MatXd& u);
    virtual VecXd hinv1_default(const MatXd& u) = 0;
    VecXd hinv2_default(const MatXd& u);

    // link between Kendall's tau and the par_bicop parameter
    VecXd tau_to_parameters(const double& tau);
    double parameters_to_tau(const VecXd& parameters);
};

#endif
