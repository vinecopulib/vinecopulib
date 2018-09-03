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
    FitControlsVinecop controls({BicopFamily::gaussian});
    controls.set_selection_criterion("loglik");

    auto u = tools_stats::simulate_uniform(20, 3);
    
    u = tools_stats::simulate_uniform(1000, 3);
    
    TVine tv(3, 3);    
    // controls.set_show_trace(true);
    tv.select_families(u, controls);
    tv.select_all(u, controls);
}

} // namespace test_tvine
