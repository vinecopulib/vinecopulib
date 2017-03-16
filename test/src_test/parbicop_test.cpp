// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/parbicop_test.hpp"

void FakeParBicopTest::set_family(BicopFamily family, int rotation)
{
    switch (family) {
        case BicopFamily::Indep:
            family_ = 0;
            break;
        case BicopFamily::Gaussian:
            family_ = 1;
            break;
        case BicopFamily::Student:
            family_ = 2;
            break;
        case BicopFamily::Clayton:
            family_ = 3;
            break;
        case BicopFamily::Gumbel:
            family_ = 4;
            break;
        case BicopFamily::Frank:
            family_ = 5;
            break;
        case BicopFamily::Joe:
            family_ = 6;
            break;
        case BicopFamily::BB1:
            family_ = 7;
            break;
        case BicopFamily::BB6:
            family_ = 8;
            break;
        case BicopFamily::BB7:
            family_ = 9;
            break;
        case BicopFamily::BB8:
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
void FakeParBicopTest::set_parameters(Eigen::VectorXd parameters)
{
    if (parameters.size() > 0)
        par_ = parameters(0);
    if (parameters.size() > 1)
        par2_ = parameters(1);
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
