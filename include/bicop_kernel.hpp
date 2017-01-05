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

#include "bicop_interpolation.hpp"
#include "bicop_class.hpp"

class KernelBicop : public Bicop {
    virtual void fit(const MatXd& data) = 0;

    VecXd pdf_default(const VecXd& u);
    VecXd hfunc1_default(const VecXd& u);
    VecXd hfunc2_default(const VecXd& u);
    VecXd hinv1_default(const VecXd& u);
    VecXd hinv2_default(const VecXd& u);

protected:
    InterpolationGrid interp_grid_;
};


#endif
