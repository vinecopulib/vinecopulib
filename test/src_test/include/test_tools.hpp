// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <fstream>

typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;
typedef Eigen::MatrixXi MatXi;

MatXd read_matxd(const char *filename, int max_buffer_size = (int) 1e6);
MatXi read_matxi(const char *filename, int max_buffer_size = (int) 1e6);
