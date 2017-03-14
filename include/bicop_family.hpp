// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <vector>

namespace families {

    enum class BicopFamily
    {
        Indep,
        Gaussian,
        Student,
        Clayton,
        Gumbel,
        Frank,
        Joe,
        BB1,
        BB6,
        BB7,
        BB8,
        TLL0
    };

    const std::vector<BicopFamily> families_parameteric = {Indep, Gaussian, Student, Clayton, Gumbel, Frank, Joe,
                                                           BB1, BB6, BB7, BB8};
    const std::vector<BicopFamily> families_nonparametric = {TLL0};
    const std::vector<BicopFamily> families_elliptical = {Gaussian, Student};
    const std::vector<BicopFamily> families_archimedean = {Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7, BB8};
    const std::vector<BicopFamily> families_rotationless = {Indep, Gaussian, Student, Frank, TLL0};
    const std::vector<BicopFamily> families_itau = {Indep, Gaussian, Student, Clayton, Gumbel, Frank, Joe};
}