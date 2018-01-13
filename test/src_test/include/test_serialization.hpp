// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "rscript.hpp"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_serialization {
using namespace vinecopulib;

TEST(serialization, bicop_serialization) {
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

TEST(serialization, vinecop_serialization) {
    // create vine with 5 variables, 2-truncated
    size_t d = 5;
    Vinecop vc(d);
    auto pc_store = Vinecop::make_pair_copula_store(d, 2);
    for (auto& tree : pc_store) {
        for (auto& pc : tree) {
            pc = Bicop(BicopFamily::bb1, 90);
        }
    }
    vc = Vinecop(pc_store, vc.get_matrix());
    
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
