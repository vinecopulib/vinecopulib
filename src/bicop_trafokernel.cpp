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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>


VecXd gaussian_kernel_1d(const VecXd& x)
{
    VecXd out(x);
    int n = x.rows();

    for (int i = 0; i < n; ++i) {
        if (std::fabs(x(i)) >= 5) {
            out(i) = 0.0;
        } else {
            out(i) = exp(- 0.5 * pow(x(i), 2)) / (sqrt(2.0 * 3.1415)) / 0.9999994267;
        }
    }
    return out;
}

VecXd gaussian_kernel_2d(const MatXd& x)
{
    return gaussian_kernel_1d(x.col(0)).cwiseProduct(gaussian_kernel_1d(x.col(1)));
}

void TrafokernelBicop::fit(const MatXd& data, __attribute__((unused)) std::string method)
{
    // default grid for evaluating the kernel estimate
    int m = 30;
    VecXd grid_points(m);
    grid_points <<
        0.000577025042390767, 0.00123962686200121, 0.00254151589573534,
        0.0049746530390053, 0.00930009771273897, 0.0166142962429879,
        0.0283788314258559, 0.0463781084506, 0.0725724690258812,
        0.108832916878933, 0.156578290481851, 0.21637844739498,
        0.287622128209223, 0.368357426310209, 0.455384362153775,
        0.544615637846225, 0.631642573689791, 0.712377871790777,
        0.78362155260502, 0.843421709518149, 0.891167083121067,
        0.927427530974119, 0.9536218915494, 0.971621168574144,
        0.983385703757012, 0.990699902287261, 0.995025346960995,
        0.997458484104265, 0.998760373137999, 0.999422974957609;


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

double TrafokernelBicop::calculate_npars()
{
    return npars_;
}
