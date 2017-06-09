// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>

namespace test_serialization {
    using namespace vinecopulib;

    TEST(serialization, bicop_serialization) {
        auto pc = Bicop(BicopFamily::bb1);
        pc.to_json("test.json");
        Bicop pc2("test.json");
        EXPECT_EQ(pc.get_rotation(), pc2.get_rotation());
        EXPECT_EQ(pc.get_family_name(), pc2.get_family_name());
        ASSERT_TRUE(pc.get_parameters().isApprox(pc2.get_parameters(), 1e-4));
    }
    
}
