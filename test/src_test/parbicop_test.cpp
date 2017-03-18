// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/parbicop_test.hpp"

void FakeParBicopTest::set_family(BicopFamily family, int rotation)
{
    switch (family) {
        case BicopFamily::indep:
            family_ = 0;
            break;
        case BicopFamily::gaussian:
            family_ = 1;
            break;
        case BicopFamily::student:
            family_ = 2;
            break;
        case BicopFamily::clayton:
            family_ = 3;
            break;
        case BicopFamily::gumbel:
            family_ = 4;
            break;
        case BicopFamily::frank:
            family_ = 5;
            break;
        case BicopFamily::joe:
            family_ = 6;
            break;
        case BicopFamily::bb1:
            family_ = 7;
            break;
        case BicopFamily::bb6:
            family_ = 8;
            break;
        case BicopFamily::bb7:
            family_ = 9;
            break;
        case BicopFamily::bb8:
            family_ = 10;
            break;
        default:
            ;
    }
    
    if (!tools_stl::is_member(family, bicop_families::rotationless)) {
        switch (rotation) {
            case 90:
                family_ += 20;
                break;
            case 180:
                family_ += 10;
                break;
            case 270:
                family_ += 30;
                break;
        }
    }
}
void FakeParBicopTest::set_parameters(std::vector<Eigen::MatrixXd> parameters)
{
    if (parameters[0].size() > 0) {
        par_ = parameters[0](0);
    }
    if (parameters[0].size() > 1) {
        par2_ = parameters[0](1);
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
