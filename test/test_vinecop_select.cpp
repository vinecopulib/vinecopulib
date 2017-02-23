/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "gtest/gtest.h"
#include "include/vinecop_class.hpp"

namespace {    
    TEST(vinecop_class, select) {
        MatXi mat(7, 7);
        mat << 4, 0, 0, 0, 0, 0, 0,
               7, 3, 0, 0, 0, 0, 0,
               3, 7, 7, 0, 0, 0, 0,
               1, 1, 5, 1, 0, 0, 0,
               2, 5, 2, 5, 2, 0, 0,
               6, 6, 1, 2, 5, 5, 0,
               5, 2, 6, 6, 6, 6, 6;
        std::vector<BicopPtr> pair_copulas(28);
        for (auto& pc : pair_copulas)
            pc = Bicop::create(3, VecXd::Constant(1, 3.0), 90);
        Vinecop vinecop(pair_copulas, mat);
        auto u = vinecop.simulate(100);  
        auto vinecop_fitted = Vinecop::structure_select(u);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
