// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace test_vinecop_discrete {

using namespace vinecopulib;

TEST(vinecop_discrete, ) {
    // simulate data and make first var discrete
    auto utmp = tools_stats::simulate_uniform(100, 10, true, {1});
    auto u = utmp;
    u.leftCols(5) = u.rightCols(5);
    u.col(0) = (utmp.col(0).array() * 2).ceil() / 2;
    u.col(5) = (utmp.col(0).array() * 2).floor() / 2;
    u.col(2) = (utmp.col(2).array() * 2).ceil() / 2;
    u.col(7) = (utmp.col(2).array() * 2).floor() / 2;

    // fit vine
    Vinecop vc(5);
    vc.set_var_types({"d", "c", "d", "c", "c"});
    FitControlsVinecop controls({BicopFamily::clayton});
    controls.set_show_trace(true);
    vc.select_all(u, controls);

    // check output
    std::cout << vc.get_rvine_structure() << std::endl; 
    auto pcs = vc.get_all_pair_copulas();
    for (size_t t = 0; t < pcs.size(); ++t) {
        for (size_t e = 0; e < pcs[t].size(); ++e) {
            std::cout << vc.get_pair_copula(t, e).str() << std::endl;
            auto vt = pcs[t][e].get_var_types();
            std::cout << vt[0] << ", " << vt[1] << std::endl;
        }
    }
}


}
