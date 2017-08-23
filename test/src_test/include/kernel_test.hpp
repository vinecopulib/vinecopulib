// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>

using namespace vinecopulib;

class TrafokernelTest : public ::testing::TestWithParam<std::string > {
public:
    TrafokernelTest();
    Bicop bicop_;
    FitControlsBicop controls;
    Eigen::MatrixXd u;
};
