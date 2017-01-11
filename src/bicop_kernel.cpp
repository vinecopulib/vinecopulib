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

KernelBicop::KernelBicop()
{
    VecXd grid_points(30);
    grid_points <<
    0.0, 0.00123962686200121, 0.00254151589573534,
    0.0049746530390053, 0.00930009771273897, 0.0166142962429879,
    0.0283788314258559, 0.0463781084506, 0.0725724690258812,
    0.108832916878933, 0.156578290481851, 0.21637844739498,
    0.287622128209223, 0.368357426310209, 0.455384362153775,
    0.544615637846225, 0.631642573689791, 0.712377871790777,
    0.78362155260502, 0.843421709518149, 0.891167083121067,
    0.927427530974119, 0.9536218915494, 0.971621168574144,
    0.983385703757012, 0.990699902287261, 0.995025346960995,
    0.997458484104265, 0.998760373137999, 1.0;
    interp_grid_ = InterpolationGrid(grid_points, MatXd::Constant(30, 30, 1.0));
}

VecXd KernelBicop::pdf_default(const MatXd& u)
{
    return interp_grid_.interpolate(u);
}
VecXd KernelBicop::hfunc1_default(const MatXd& u)
{
    return interp_grid_.intergrate_1d(u, 1);
}
VecXd KernelBicop::hfunc2_default(const MatXd& u)
{
    return interp_grid_.intergrate_1d(u, 2);
}
VecXd KernelBicop::hinv1_default(const MatXd& u)
{
    return hinv1_num(u);
}
VecXd KernelBicop::hinv2_default(const MatXd& u)
{
    return hinv2_num(u);
}

// TODO
double KernelBicop::parameters_to_tau(const VecXd __attribute__((unused))& parameters)
{
    throw std::runtime_error(
        "parameters_to_tau not yet implemented for kernel estimator"
    );
}

double KernelBicop::calculate_npars()
{
    return npars_;
}
