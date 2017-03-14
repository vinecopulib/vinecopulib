// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>

typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;
typedef Eigen::MatrixXi MatXi;
typedef Eigen::VectorXi VecXi;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatXb;

VecXd invert_f(
    const VecXd &x, std::function<VecXd(const VecXd&)> f,
    const double lb = 1e-20,
    const double ub = 1-1e-20,
    int n_iter = 35
);
