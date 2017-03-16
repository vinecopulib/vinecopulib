// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "tools_bicopselect.hpp"
#include "tools_stats.hpp"
#include "tools_stl.hpp"
#include <cmath>

std::vector<double> get_c1c2(const Eigen::MatrixXd& data, double tau)
{
    int n = data.rows();
    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, 2);
    Eigen::MatrixXd z1 = x;
    Eigen::MatrixXd z2 = x;
    x = tools_stats::qnorm(data);
    
    int count1 = 0, count2 = 0;
    for (int j = 0; j < n; ++j) {
        if (tau > 0) {
            if ((x(j, 0) > 0) && (x(j, 1) > 0)) {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if ((x(j, 0) < 0) && (x(j, 1) < 0)) {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        } else {
            if ((x(j, 0) < 0) && (x(j, 1) > 0)) {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if ((x(j, 0) > 0) && (x(j, 1) < 0)) {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        }
    }

    return  {
        tools_stats::pairwise_cor(z1.block(0, 0, count1 - 1, 2)),
        tools_stats::pairwise_cor(z2.block(0, 0, count2 - 1, 2))
    };
}

bool preselect_family(
    std::vector<double> c, 
    double tau, 
    vinecopulib::BicopFamily family, 
    int rotation, 
    bool is_rotationless
)
{
    using namespace vinecopulib;
    using namespace tools_stl;
    int c1 = c[0];
    int c2 = c[1];
    bool preselect = false;
    if (is_rotationless) {
        preselect = true;
        if ((std::fabs(c1 - c2) > 0.3) & (family == BicopFamily::Frank))
            preselect = false;
    } else {
        if (is_member(family, bicop_families::BB)) {
            if ((tau > 0) && is_member(rotation, {0, 180})) {
                preselect = true;
            }
            if ((tau < 0) && is_member(rotation, {90, 270})) {
                preselect = true;
            }
        }
        bool is_90or180 = is_member(rotation, {90, 180});
        if (c1 - c2 > 0.05) {
            if (is_member(family, bicop_families::lt) & is_90or180) { 
                preselect = true;
            }
            if (is_member(family, bicop_families::ut) & !is_90or180) {
                preselect = true;
            }
        } else if (c1 - c2 < -0.05) {
            if (is_member(family, bicop_families::lt) & !is_90or180) { 
                preselect = true;
            }
            if (is_member(family, bicop_families::ut) & is_90or180) {
                preselect = true;
            }
        } else {
            if ((tau > 0) && is_member(rotation, {0, 180})) {
                preselect = true;
            }
            if ((tau < 0) && is_member(rotation, {90, 270})) {
                preselect = true;
            }
        }
    }
    
    return preselect;
}
