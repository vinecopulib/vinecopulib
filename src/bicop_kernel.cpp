/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "bicop_kernel.hpp"

KernelBicop::KernelBicop()
{
    // construct default grid (equally spaced on Gaussian scale)
    int m = 30;
    VecXd grid_points(m);
    for (int i = 0; i < m; ++i)
        grid_points(i) = - 3.25 + i * (6.25 / (double) m);
    grid_points = pnorm(grid_points);

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
double KernelBicop::parameters_to_tau(const VecXd &)
{
    throw std::runtime_error(
        "parameters_to_tau not yet implemented for kernel estimator"
    );
}

double KernelBicop::calculate_npars()
{
    return npars_;
}

void KernelBicop::flip()
{
    interp_grid_.flip();
}
