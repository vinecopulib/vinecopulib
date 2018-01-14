// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_eigen.hpp>
#include "rscript.hpp"

class VinecopTest : public ::testing::Test {
public:
    VinecopTest();

    Eigen::MatrixXd u;
    Eigen::VectorXd f;
    Eigen::MatrixXd sim;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> model_matrix;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> vc_matrix;
};
