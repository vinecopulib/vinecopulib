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
    // hfunction and its inverse
    VecXd hfunc2(const MatXd& u);
    VecXd hinv2(const MatXd& u);

    // link between Kendall's tau and the par_bicop parameter
    VecXd tau_to_par(const double& tau);
    double par_to_tau(const VecXd& parameters);
};

#endif
