// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "rscript.hpp"
#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

namespace test_serialization {
using namespace vinecopulib;

TEST(serialization, bicop_serialization)
{
  auto pc = Bicop(BicopFamily::bb1);
  pc.to_json("temp");
  Bicop pc2("temp");

  // Remove temp file
  std::string cmd = rm + "temp";
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  EXPECT_EQ(pc.get_rotation(), pc2.get_rotation());
  EXPECT_EQ(pc.get_family_name(), pc2.get_family_name());
  ASSERT_TRUE(pc.get_parameters().isApprox(pc2.get_parameters(), 1e-4));
}

TEST(serialization, vinecop_serialization)
{

  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  // create vine with 7 variables, 2-truncated
  size_t d = 7;
  auto pc_store = Vinecop::make_pair_copula_store(d, 5);
  for (auto& tree : pc_store) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::bb1, 90);
    }
  }

  auto vc = Vinecop(pc_store, mat);

  // serialize
  vc.to_json("temp");

  // unserialize
  auto vc2 = Vinecop("temp");

  // Remove temp file
  std::string cmd = rm + "temp";
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  EXPECT_EQ(vc.get_all_rotations(), vc2.get_all_rotations());
  EXPECT_EQ(vc.get_all_families(), vc2.get_all_families());
}
}
