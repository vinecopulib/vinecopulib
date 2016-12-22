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

#ifndef VINECOPLIB_BICOP_CLAYTON_HPP
#define VINECOPLIB_BICOP_CLAYTON_HPP

#include "bicop_archimedean.hpp"

class ClaytonBicop : public ArchimedeanBicop {

public:
    // constructor
    ClaytonBicop();
    ClaytonBicop(const VecXd& parameters);
    ClaytonBicop(const VecXd& parameters, const int& rotation);

    // generator, its inverse and derivatives for the archimedean copula
    VecXd generator(const VecXd& u);
    VecXd generator_inv(const VecXd& u);
    VecXd generator_derivative(const VecXd& u);
    VecXd generator_derivative2(const VecXd& u);

    // bounds on the copula parameter
    MatXd get_bounds_standard();

    // inverse hfunction
    VecXd hinv(const MatXd& u);

    // link between Kendall's tau and the par_bicop parameter
    VecXd tau_to_par(const double& tau);
    double par_to_tau(const VecXd& parameters);
};


#endif
