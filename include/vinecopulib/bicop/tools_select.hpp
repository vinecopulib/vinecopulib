// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vector>

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/bicop/class.hpp>

namespace vinecopulib {
namespace tools_select {

std::vector <Bicop> create_candidate_bicops(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
    const FitControlsBicop &controls);

std::vector <BicopFamily> get_candidate_families(
    const FitControlsBicop &controls);

void preselect_candidates(
    std::vector <Bicop> &bicops,
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
    double tau);

std::vector<double> get_c1c2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &data, double tau);

bool preselect_family(std::vector<double> c, double tau, const Bicop &bicop);

}
}

#include <vinecopulib/bicop/implementation/tools_select.ipp>
