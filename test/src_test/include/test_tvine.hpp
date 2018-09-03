// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "test_tvine.hpp"
#include <vinecopulib/vinecop/tvine.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace test_tvine
{

using namespace vinecopulib;

TEST(test_tvine, playground)
{
    auto u = tools_stats::simulate_uniform(10, 3);
    Vinecop vc(u);
    TVine tv(vc.get_rvine_structure(), 3);
    u = tools_stats::simulate_uniform(10, 9);
    tv.select_families(u);
}

} // namespace test_tvine
