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
#include <include/boost_tools.hpp>

int main(int argc, char **argv) {

    VecXd x = VecXd::Zero(3);
    x(0) = -0.5;
    x(2) = 0.5;
    std::cout << dnorm(x) << std::endl;
    std::cout << pnorm(x) << std::endl;
    x(0) = 0.05;
    x(1) = 0.5;
    x(2) = 0.95;
    std::cout << qnorm(x) << std::endl;

    return 0;
}
