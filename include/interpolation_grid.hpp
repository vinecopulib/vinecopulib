// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "tools_eigen.hpp"

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
