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

#include "bicop_trafokernel.hpp"
#include "bicop.hpp"

TrafokernelBicop::TrafokernelBicop()
{
    family_ = 1001;
    family_name_ = "Transformation kernel";
    rotation_ = 0;
    association_direction_ = "both";
}

VecXd gaussian_kernel_2d(const MatXd& x)
{
    return dnorm(x).rowwise().prod();
}

void TrafokernelBicop::fit(const MatXd& data, std::string)
{
    // construct default grid (equally spaced on Gaussian scale)
    int m = 30;
    VecXd grid_points(m);
    for (int i = 0; i < m; ++i)
        grid_points(i) = - 3.25 + i * (6.25 / (double) m);
    grid_points = pnorm(grid_points);

    // expand the interpolation grid; a matrix with two columns where each row
    // contains one combination of the grid points
    MatXd grid_2d(m * m, 2);
    int k = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            grid_2d(k, 0) = grid_points(i);
            grid_2d(k, 1) = grid_points(j);
            ++k;
        }
    }

    // transform evaluation grid and data by inverse Gaussian cdf
    MatXd z = qnorm(grid_2d);
    MatXd z_data = qnorm(data);

    // apply normal density to z (used later for normalization)
    MatXd phi = dnorm(z);
    MatXd phi_data = dnorm(z_data);

    // find bandwidth matrix
    int n = data.rows();
    MatXd centered = z_data.rowwise() - z_data.colwise().mean();
    MatXd cov = (centered.adjoint() * centered) / double(n - 1);

    Eigen::SelfAdjointEigenSolver<MatXd> takes_root(cov);
    MatXd cov_root = takes_root.operatorSqrt();
    MatXd B =  1.25 * std::pow(n, - 1.0 / 6.0) * cov_root.transpose();

    // apply bandwidth matrix
    z = (B.inverse() * z.transpose()).transpose();
    z_data = (B.inverse() * z_data.transpose()).transpose();

    // compute estimator on each evaluation point
    VecXd kernels(n);
    double det_B = B.determinant();
    int i = 0;
    int j = 0;
    MatXd values(m, m);
    for (int k = 0; k < m * m; ++k) {
        kernels = gaussian_kernel_2d((z_data - z.row(k).replicate(n, 1)));
        values(i, j) = kernels.mean() / (det_B * phi.row(k).prod());
        ++i;
        if (i % m == 0) {
            i = 0;
            ++j;
        }
    }

    // for interpolation, we shift the limiting gridpoints to 0 and 1
    grid_points(0) = 0.0;
    grid_points(m - 1) = 1.0;
    interp_grid_ = InterpolationGrid(grid_points, values);
    
    // compute effective number of parameters
    double K0 = gaussian_kernel_2d(MatXd::Constant(1, 2, 0.0))(0);
    VecXd scale = phi_data.rowwise().prod();
    npars_ =  K0 / det_B / (scale.array() * this->pdf(data).array()).mean();
}
