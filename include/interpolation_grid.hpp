// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "tools_eigen.hpp"

namespace vinecopulib
{
    //! A class for cubic spline interpolation of bivariate copulas
    //!
    //! The class is used for implementing kernel estimators. It makes storing the
    //! observations obsolete and allows for fast numerical integration.
    class InterpolationGrid {
    public:
        InterpolationGrid() {}
        InterpolationGrid(const VectorXd& grid_points, const MatrixXd& values);

        void flip();

        VectorXd interpolate(const MatrixXd& x);
        VectorXd intergrate_1d(const MatrixXd& u, const int& cond_var);
        VectorXd inv_intergrate_1d(const MatrixXd& u, const int& cond_var);

    private:
        // Utility functions for spline Interpolation
        double cubic_poly(const double& x, const VectorXd& a);
        double cubic_indef_integral(const double& x, const VectorXd& a);
        double cubic_integral(const double& lower, const double& upper, const VectorXd& a);
        double inv_cubic_integral(const double& q, const VectorXd& a);
        VectorXd find_coefs(const VectorXd& vals, const VectorXd& grid);
        double interp_on_grid(const double& x, const VectorXd& vals, const VectorXd& grid);

        // Utility functions for integration
        double int_on_grid(const double& upr, const VectorXd& vals, const VectorXd& grid);
        double inv_int_on_grid(const double& qq, const VectorXd& vals, const VectorXd& grid);

        VectorXd grid_points_;
        MatrixXd values_;
    };
}
