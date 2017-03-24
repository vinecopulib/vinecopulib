// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/trafokernel_test.hpp"
#include "include/misc/tools_stats.hpp"

TrafokernelTest::TrafokernelTest() {
    bicop_ = Bicop(vinecopulib::BicopFamily::tll0, 0);
    u = tools_stats::simulate_uniform(100, 2);
}
