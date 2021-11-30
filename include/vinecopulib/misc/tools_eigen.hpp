// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <functional>
#include <vinecopulib/misc/tools_linalg.hpp>

namespace vinecopulib {

//! Tools for working with Eigen types
namespace tools_eigen {
//! An `Eigen::Matrix` containing `bool`s (similar to `Matrix`).
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

template<typename T>
Matrix
unaryExpr_or_nan(const Matrix& x, const T& func)
{
  return x.unaryExpr([&func](const double& y) {
    if ((boost::math::isnan)(y)) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      return func(y);
    }
  });
}

template<typename T>
Vector
binaryExpr_or_nan(const Matrix& u, const T& func)
{
  auto func_or_nan = [&func](const double& u1, const double& u2) {
    if ((boost::math::isnan)(u1) | (boost::math::isnan)(u2)) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      return func(u1, u2);
    }
  };
  return u.col(0).binaryExpr(u.col(1), func_or_nan);
}

void
remove_nans(Matrix& x);

void
remove_nans(Matrix& x, Vector& weights);

void
trim(Matrix& x, const double& lower = 1e-10, const double& upper = 1 - 1e-10);

void
trim(Vector& x,
     const double& lower = 1e-10,
     const double& upper = 1 - 1e-10);

bool
check_if_in_unit_cube(const Matrix& u);

Matrix
swap_cols(Matrix u);

Vector
unique(const Vector& x);

Vector
invert_f(const Vector& x,
         std::function<Vector(const Vector&)> f,
         const double lb = 1e-20,
         const double ub = 1 - 1e-20,
         int n_iter = 35);

Matrix
expand_grid(const Vector& grid_points);

Matrix
read_matxd(const char* filename, int max_buffer_size = static_cast<int>(1e6));

Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
read_matxs(const char* filename, int max_buffer_size = static_cast<int>(1e6));
}
}

#include <vinecopulib/misc/implementation/tools_eigen.ipp>
