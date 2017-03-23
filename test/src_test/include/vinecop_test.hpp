// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "include/misc/tools_eigen.hpp"
#include "rscript.hpp"

using namespace vinecopulib;

class VinecopTest : public ::testing::Test {
public:
    VinecopTest();
    Eigen::MatrixXd u;
    Eigen::VectorXd f;
    Eigen::MatrixXd sim;
    Eigen::MatrixXi model_matrix;
    Eigen::MatrixXi vc_matrix;
};
