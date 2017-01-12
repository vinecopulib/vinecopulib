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

#include "bicop_trafokernel.hpp"
#include "bicop.hpp"

TrafokernelBicop::TrafokernelBicop()
{
    family_ = 1001;
    rotation_ = 0;
    association_direction_ = "both";
}

VecXd gaussian_kernel_2d(const MatXd& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_ran_ugaussian_pdf)).rowwise().prod();
}

void TrafokernelBicop::fit(const MatXd& data, __attribute__((unused)) std::string method)
{
    // construct default grid (equally spaced on Gaussian scale)
    int m = 30;
    VecXd grid_points(m);
    for (int i = 0; i < m; ++i)
        grid_points(i) = gsl_cdf_ugaussian_P(- 3.25 + i * (6.25 / (double) m));

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
    MatXd z = grid_2d.unaryExpr(std::ptr_fun(gsl_cdf_ugaussian_Pinv));
    MatXd z_data = data.unaryExpr(std::ptr_fun(gsl_cdf_ugaussian_Pinv));

    // apply normal density to z (used later for normalization)
    MatXd phi = z.unaryExpr(std::ptr_fun(gsl_ran_ugaussian_pdf));
    MatXd phi_data = z_data.unaryExpr(std::ptr_fun(gsl_ran_ugaussian_pdf));

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
