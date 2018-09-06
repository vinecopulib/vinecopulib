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

    auto u = tools_stats::simulate_uniform(30, 4);
        
    TVine tv(4, 0);
    
    auto cs_struct = tv.get_cs_structure();
    // // cs_struct.truncate(2);

    TVine tv2(cs_struct, 1, 2, 4);
    
    controls.set_show_trace(true);
    // tv2.select_families(u, controls);
    tv.select_all(u, controls);
    // 
    // TVine(tv.get_all_pair_copulas(), tv.get_rvine_structure(), 2);
    // std::cout << tv.simulate(10) << std::endl;

}

} // namespace test_tvine
