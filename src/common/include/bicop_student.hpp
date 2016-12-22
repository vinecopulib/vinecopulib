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

#ifndef VINECOPLIB_STUDENT_BICOP_HPP
#define VINECOPLIB_STUDENT_BICOP_HPP

#include "bicop_elliptical.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

class StudentBicop : public EllipticalBicop {

public:
    // constructor
    StudentBicop();
    StudentBicop(const VecXd& parameters);
    StudentBicop(const VecXd& parameters, const int& rotation);

    // PDF
    VecXd pdf(const MatXd& u);

    // hfunction
    VecXd hfunc1(const MatXd& u);

    // inverse hfunction
    VecXd hinv1(const MatXd& u);
};


#endif
