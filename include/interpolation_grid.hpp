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

#ifndef VINECOPULIB_BICOP_INTERPOLATION_HPP
#define VINECOPULIB_BICOP_INTERPOLATION_HPP

#include <Eigen/Dense>

typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

//! A class for cubic spline interpolation of bivariate copulas
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration.
class InterpolationGrid {
public:
    InterpolationGrid() {}
    InterpolationGrid(const VecXd& grid_points, const MatXd& values);

    void flip();

    VecXd interpolate(const MatXd& x);
    VecXd intergrate_1d(const MatXd& u, const int& cond_var);
    VecXd inv_intergrate_1d(const MatXd& u, const int& cond_var);

private:
    // Utility functions for spline Interpolation
    double cubic_poly(const double& x, const VecXd& a);
    double cubic_indef_integral(const double& x, const VecXd& a);
    double cubic_integral(const double& lower, const double& upper, const VecXd& a);
    double inv_cubic_integral(const double& q, const VecXd& a);
    VecXd find_coefs(const VecXd& vals, const VecXd& grid);
    double interp_on_grid(const double& x, const VecXd& vals, const VecXd& grid);

    // Utility functions for integration
    double int_on_grid(const double& upr, const VecXd& vals, const VecXd& grid);
    double inv_int_on_grid(const double& qq, const VecXd& vals, const VecXd& grid);

    VecXd grid_points_;
    MatXd values_;
};

#endif
