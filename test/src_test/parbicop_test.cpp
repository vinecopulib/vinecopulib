// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/parbicop_test.hpp"

void FakeParBicopTest::set_family(int family, int rotation)
{
    std::vector<int> rotation_less_fams = {0, 1, 2, 5};
    if (tools_stl::is_member(family, rotation_less_fams) || rotation == 0)
    {
        family_ = family;
    }
    else
    {
        if (rotation == 90)
        {
            family_ = family + 20;
        }
        else if (rotation == 180)
        {
            family_ = family + 10;
        }
        else
        {
            family_ = family + 30;
        }
    }
}
void FakeParBicopTest::set_parameters(Eigen::VectorXd parameters)
{
    if (parameters.size() > 0)
    {
        par_ = parameters(0);
    }
    if (parameters.size() > 1)
    {
        par2_ = parameters(1);
    }
}
void FakeParBicopTest::set_n(int n)
{
    n_ = n;
}
int FakeParBicopTest::get_family()
{
    return family_;
}
int FakeParBicopTest::get_n()
{
    return n_;
}
double FakeParBicopTest::get_par()
{
    return par_;
}
double FakeParBicopTest::get_par2()
{
    return par2_;
}
