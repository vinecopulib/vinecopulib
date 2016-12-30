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
#include "include/bicop.hpp"
#include <vector>

namespace {
    // Test if the C++ implementation of the par_to_tau and tau_to_par is correct
    TEST(ParBicopTest, bicop_creates_right_copula) {
        std::vector<int> fams = {0, 1, 2, 3, 4, 5, 6};
        std::vector<int> rots = {0, 90, 180, 270};

        for (unsigned int i = 0; i < fams.size(); ++i) {
            for (unsigned int j = 0; j < rots.size(); ++j) {
                Bicop_ptr bicop = Bicop::create(fams[i], VecXd::Zero(1), rots[j]);
                EXPECT_EQ(bicop->get_family(), fams[i]);
                EXPECT_EQ(bicop->get_rotation(), rots[j]);
                EXPECT_EQ(bicop->get_parameters(), VecXd::Zero(1));
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
