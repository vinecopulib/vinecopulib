// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/vinecop/tvine.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace test_tvine
{

using namespace vinecopulib;

TEST(test_tvine, playground)
{
    FitControlsVinecop controls({BicopFamily::clayton});
    
    auto u = tools_stats::simulate_uniform(30, 4);
        
    TVine tv(4, 0);
    
    RVineStructure cs_struct = tv.get_cs_structure();
    // cs_struct.truncate(2);

    
    TVine tv2(cs_struct, 1, 1, 3);
    std::cout << cs_struct << std::endl;
    std::cout << tv2.get_rvine_structure() << std::endl;
    // 
    controls.set_show_trace(true);
    tv2.select_all(u, controls);
    std::cout << tv2.get_tvine_structure() << std::endl;

    // std::cout << tv2.simulate(10) << std::endl;
    // 
    // // tv2.select_all(u, controls);
    // // std::cout << tv2.get_rvine_structure() << std::endl;

    TVine tv3(tv2.get_all_pair_copulas(), 
              tv2.get_cs_structure(),
              tv2.get_p(),
              tv2.get_in_vertex(),
              tv2.get_out_vertex());
    std::cout << tv3.simulate(10) << std::endl;
}

} // namespace test_tvine
