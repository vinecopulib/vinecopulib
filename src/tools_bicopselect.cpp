// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "tools_bicopselect.hpp"
#include "tools_stats.hpp"

std::vector<double> get_c1c2(const MatXd& data, double tau)
{
    int n = data.rows();
    MatXd x = MatXd::Zero(n,2);
    MatXd z1 = x;
    MatXd z2 = x;
    x = qnorm(data);
    int count1 = 0, count2 = 0;
    for (int j = 0; j < n; ++j) {
        if (tau > 0)
        {
            if (x(j,0) > 0 && x(j,1) > 0)
            {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if (x(j,0) < 0 && x(j,1) < 0)
            {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        } else {
            if (x(j,0) < 0 && x(j,1) > 0)
            {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if (x(j,0) > 0 && x(j,1) < 0)
            {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        }
    }

    std::vector<double> c = {
            pairwise_cor(z1.block(0,0,count1-1,2)),
            pairwise_cor(z2.block(0,0,count2-1,2))
    };
    return c;
}

bool preselect_family(double c1, double c2, double tau, int family, int rotation, bool is_rotationless)
{
    bool preselect = false;
    if (is_rotationless)
    {
        if (std::fabs(c1-c2) > 0.3)
        {
            if (family == 2 || family == 0 || family == 1001)
            {
                preselect = true;
            }
        } else
        {
            preselect = true;
        }
    } else
    {
        bool is_90or180 = (rotation == 90 || rotation == 180);
        if (family == 7 || family == 8 || family == 9 || family == 10)
        {
            if (tau > 0 && (rotation == 0 || rotation == 180))
            {
                preselect = true;
            }
            if (tau < 0 && (rotation == 90 || rotation == 270))
            {
                preselect = true;
            }
        }
        if (c1 - c2 > 0.05)
        {
            if (family == 3)
            {
                if (is_90or180)
                {
                    preselect = true;
                }
            }
            if (family == 4 || family == 6)
            {
                if (!is_90or180)
                {
                    preselect = true;
                }
            }
        } else if (c1 - c2 < -0.05)
        {
            if (family == 3)
            {
                if (!is_90or180)
                {
                    preselect = true;
                }
            }
            if (family == 4 || family == 6)
            {
                if (is_90or180)
                {
                    preselect = true;
                }
            }
        } else {
            if (tau > 0 && (rotation == 0 || rotation == 180))
            {
                preselect = true;
            }
            if (tau < 0 && (rotation == 90 || rotation == 270))
            {
                preselect = true;
            }
        }
    }
    return preselect;
}