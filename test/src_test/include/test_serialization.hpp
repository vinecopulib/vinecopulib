// Copyright © 2017 Thomas Nagler and Thibault Vatter
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

        size_t d = 5;
        auto vc = Vinecop(d);
        vc.to_json("temp");
        auto vc2 = Vinecop("temp");

        // Remove temp file
        std::string cmd = rm + "temp";
        int sys_exit_code = system(cmd.c_str());
        if (sys_exit_code != 0) {
            throw std::runtime_error("error in system call");
        }

        for (size_t tree = 0; tree < d - 1; ++tree) {
            for (size_t edge = 0; edge < d - tree - 1; ++edge) {
                EXPECT_EQ(vc2.get_rotation(tree, edge), 0);
                EXPECT_EQ(vc2.get_family(tree, edge), BicopFamily::indep);
            }
        }
    }
    
}
