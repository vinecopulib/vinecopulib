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

#ifndef VINECOPULIB_BICOP_KERNEL_HPP
#define VINECOPULIB_BICOP_KERNEL_HPP

#include "interpolation_grid.hpp"
#include "bicop_class.hpp"

class KernelBicop : public Bicop {
public:
    KernelBicop();

    VecXd pdf_default(const MatXd& u);
    VecXd hfunc1_default(const MatXd& u);
    VecXd hfunc2_default(const MatXd& u);
    VecXd hinv1_default(const MatXd& u);
    VecXd hinv2_default(const MatXd& u);

    double parameters_to_tau(const VecXd __attribute__((unused))& parameters);
    double calculate_npars();

protected:
    InterpolationGrid interp_grid_;
    double npars_;
};


#endif
