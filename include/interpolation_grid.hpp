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

#pragma once

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

