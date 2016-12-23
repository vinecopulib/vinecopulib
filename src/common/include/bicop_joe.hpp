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

#ifndef VINECOPULIB_BICOP_JOE_HPP
#define VINECOPULIB_BICOP_JOE_HPP


#include "bicop_archimedean.hpp"
#include "gsl/gsl_sf_psi.h"

class JoeBicop : public ArchimedeanBicop {

public:
    // constructor
    JoeBicop();
    JoeBicop(const VecXd& parameters);
    JoeBicop(const VecXd& parameters, const int& rotation);

    // generator, its inverse and derivatives for the archimedean copula
    VecXd generator(const VecXd& u);
    VecXd generator_inv(const VecXd& u);
    VecXd generator_derivative(const VecXd& u);
    VecXd generator_derivative2(const VecXd& u);

    // inverse hfunction
    VecXd hinv(const MatXd& u);

    // link between Kendall's tau and the par_bicop parameter
    VecXd tau_to_par(const double& tau);
    double par_to_tau(const VecXd& par);
};

double qcondjoe(double* q, double* u, double* de);


#endif
