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
#include "include/rvine_matrix.hpp"

namespace {
    TEST(rvine_matrix, can_convert_to_natural_order) {
        MatXi mat(5, 5);
        mat << 5, 0, 0, 0, 0,
               2, 2, 0, 0, 0,
               3, 3, 3, 0, 0,
               1, 4, 4, 4, 0,
               4, 1, 1, 1, 1;
        MatXi no_mat(5, 5);
        no_mat << 5, 0, 0, 0, 0,
                  4, 4, 0, 0, 0,
                  3, 3, 3, 0, 0,
                  1, 2, 2, 2, 0,
                  2, 1, 1, 1, 1;
        MatXi my_no = RVineMatrix::to_natural_order(mat);
        EXPECT_EQ(my_no, no_mat);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
