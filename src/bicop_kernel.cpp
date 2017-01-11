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
    // construct default grid (equally spaced on Gaussian scale)
    int m = 30;
    VecXd grid_points(m);
    for (int i = 0; i < m; ++i)
        grid_points(i) = gsl_cdf_ugaussian_P(- 3.25 + i * (6.25 / (double) m));

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
