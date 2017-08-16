// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_os.hpp>
#include <stan/command/stanc_helper.hpp>

using namespace vinecopulib;

namespace test_stan {

    TEST(test_stan, stanc_works) {
        std::stringstream out;
        std::stringstream err;
        int argc = 4;
        std::vector<const char*> argv_vec;
        argv_vec.push_back("main");
        argv_vec.push_back("--name=rosenbrock");
        argv_vec.push_back("--o=bin/rosenbrock.hpp");
        argv_vec.push_back("bin/rosenbrock.stan");
        const char** argv = &argv_vec[0];
        int rc = stanc_helper(argc, argv, &out, &err);
        EXPECT_TRUE(rc == 0) << "error=" << err.str() << std::endl;
    }

}
