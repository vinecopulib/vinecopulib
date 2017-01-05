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

#include "bicop_kernel.hpp"

VecXd KernelBicop::pdf_default(const VecXd& u)
{
    return interp_grid_.interpolate(u);
}
VecXd KernelBicop::hfunc1_default(const VecXd& u)
{
    return interp_grid_.intergrate_1d(u, 1);
}
VecXd KernelBicop::hfunc2_default(const VecXd& u)
{
    return interp_grid_.intergrate_1d(u, 2);
}
VecXd KernelBicop::hinv1_default(const VecXd& u)
{
    return interp_grid_.inv_intergrate_1d(u, 1);
}
VecXd KernelBicop::hinv2_default(const VecXd& u)
{
    return interp_grid_.inv_intergrate_1d(u, 2);
}
