// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_linalg.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! A class for cubic spline interpolation of bivariate copulas
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration.
class InterpolationGrid
{
public:
  InterpolationGrid() {}

  InterpolationGrid(const Vector& grid_points,
                    const Matrix& values,
                    int norm_times = 3);

  Matrix get_values() const;

  void set_values(const Matrix& values, int norm_times = 3);

  void flip();

  void normalize_margins(int times);

  Vector interpolate(const Matrix& x);

  Vector integrate_1d(const Matrix& u, size_t cond_var);

  Vector integrate_2d(const Matrix& u);

private:
  Eigen::Matrix<ptrdiff_t, 1, 2> get_indices(double x0, double x1);
  double bilinear_interpolation(double z11,
                                double z12,
                                double z21,
                                double z22,
                                double x1,
                                double x2,
                                double y1,
                                double y2,
                                double x,
                                double y);
  double int_on_grid(const double& upr,
                     const Vector& vals,
                     const Vector& grid);

  Vector grid_points_;
  Matrix values_;
};
}
}

#include <vinecopulib/misc/implementation/tools_interpolation.ipp>
