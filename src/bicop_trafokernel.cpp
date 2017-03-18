// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_trafokernel.hpp"
#include "bicop.hpp"
#include "tools_stats.hpp"

namespace vinecopulib
{
    TrafokernelBicop::TrafokernelBicop()
    {
        family_ = BicopFamily::tll0;
        rotation_ = 0;
        parameters_ = {
            Eigen::VectorXd::Constant(1, 1),  // multiplier
            Eigen::Matrix2d::Identity()       // bandwidth matrix
        };
    }

    Eigen::VectorXd gaussian_kernel_2d(const Eigen::MatrixXd& x)
    {
        return tools_stats::dnorm(x).rowwise().prod();
    }

    void TrafokernelBicop::fit(const Eigen::MatrixXd& data, std::string)
    {
        // construct default grid (equally spaced on Gaussian scale)
        int m = 30;
        Eigen::VectorXd grid_points(m);
        for (int i = 0; i < m; ++i)
            grid_points(i) = - 3.25 + i * (6.25 / (double) m);
        grid_points = tools_stats::pnorm(grid_points);

        // expand the interpolation grid; a matrix with two columns where each row
        // contains one combination of the grid points
        Eigen::MatrixXd grid_2d(m * m, 2);
        int k = 0;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                grid_2d(k, 0) = grid_points(i);
                grid_2d(k, 1) = grid_points(j);
                ++k;
            }
        }

        // transform evaluation grid and data by inverse Gaussian cdf
        Eigen::MatrixXd z = tools_stats::qnorm(grid_2d);
        Eigen::MatrixXd z_data = tools_stats::qnorm(data);

        // apply normal density to z (used later for normalization)
        Eigen::MatrixXd phi = tools_stats::dnorm(z);
        Eigen::MatrixXd phi_data = tools_stats::dnorm(z_data);

        // find bandwitools_stats::dth matrix
        int n = data.rows();
        Eigen::MatrixXd centered = z_data.rowwise() - z_data.colwise().mean();
        Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(n - 1);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> takes_root(cov);
        Eigen::MatrixXd cov_root = takes_root.operatorSqrt();
        Eigen::MatrixXd B =  1.25 * std::pow(n, - 1.0 / 6.0) * cov_root.transpose();
        parameters_[1] = B;

        // apply bandwitools_stats::dth matrix
        z = (B.inverse() * z.transpose()).transpose();
        z_data = (B.inverse() * z_data.transpose()).transpose();

        // compute estimator on each evaluation point
        Eigen::VectorXd kernels(n);
        double det_B = B.determinant();
        int i = 0;
        int j = 0;
        Eigen::MatrixXd values(m, m);
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
        double K0 = gaussian_kernel_2d(Eigen::MatrixXd::Constant(1, 2, 0.0))(0);
        Eigen::VectorXd scale = phi_data.rowwise().prod();
        npars_ =  K0 / det_B / (scale.array() * this->pdf(data).array()).mean();
    }
}
