// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "r_parity.hpp"
#include "gtest/gtest.h"
#include <vinecopulib.hpp>

//! @brief Fixture whose data comes from VineCopula, via the R oracle.
//!
//! Every test on this fixture consumes that data, so the whole fixture skips
//! when the build has no Rscript; loading happens in SetUp() because
//! GTEST_SKIP() is not usable from a constructor.
class VinecopTest : public ::testing::Test
{
public:
  void SetUp() override;

  Eigen::MatrixXd u;
  Eigen::VectorXd f;
  Eigen::MatrixXd sim;
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> model_matrix;
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> vc_matrix;
};
