/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecoplib.

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

#ifndef VINECOPLIB_NORMAL_BICOP_H
#define VINECOPLIB_NORMAL_BICOP_H

#include "bicop_elliptical.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
//#include <boost/math/distributions/normal.hpp>

class NormalBicop : public EllipticalBicop {

public:
    // constructor
    NormalBicop();
    NormalBicop(double rho);

    // bounds on the copula parameter
    MatXd get_bounds();

    // hfunction
    VecXd hfunc1(const MatXd &u);

    // inverse hfunction
    VecXd hinv1(const MatXd &u);

    // PDF
    VecXd pdf(const MatXd &u);
};


#endif //VINECOPLIB_NORMAL_BICOP_H
