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

#ifndef VINECOPLIB_BICOP_H_INCLUDED__
#define VINECOPLIB_BICOP_H_INCLUDED__

#include <Eigen/Dense>
#include <random>
#include <nlopt.hpp>
#include <iostream>
typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

extern "C" {
#include "ktau.h"
}

class Bicop
{

public:

    // hfunctions: the conditioning variable is put second
    virtual VecXd hfunc1(const MatXd &u) = 0;
    virtual VecXd hfunc2(const MatXd &u) = 0;

    // PDF
    virtual VecXd pdf(const MatXd &u) = 0;

    // Inverse of the h-functions
    virtual VecXd hinv1(const MatXd &u) = 0;
    virtual VecXd hinv2(const MatXd &u) = 0;
    VecXd hinv1_num(const MatXd &u);
    VecXd hinv2_num(const MatXd &u);

    // Number of parameters
    virtual int calculate_npars() = 0;

    // Simulation
    MatXd simulate(int n);

    // fit statistics
    double loglik(MatXd &u);
    double aic(MatXd &u);
    double bic(MatXd &u);

};

#endif // __BICOP_H_INCLUDED__
