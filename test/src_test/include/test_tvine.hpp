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

    TVine tv2(cs_struct, 1, 4, 2);

    tv2.select_families(u, controls);
    std::cout << tv2.get_rvine_structure() << std::endl;

    std::cout << tv2.simulate(10) << std::endl;
    std::cout << tv2.simulate_conditional(10, u) << std::endl;
    std::cout << tv2.simulate_ahead(10, u) << std::endl;

    tv2.select_all(u, controls);
    std::cout << tv2.get_rvine_structure() << std::endl;

}

} // namespace test_tvine
