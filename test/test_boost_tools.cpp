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
#include "include/boost_tools.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <Eigen/Dense>
#include <ctime>

template<typename T> T dnorm_gsl(const T& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_ran_ugaussian_pdf));
};
template<typename T> T pnorm_gsl(const T& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_cdf_ugaussian_P));
};
template<typename T> T qnorm_gsl(const T& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_cdf_ugaussian_Pinv));
};

int main(int argc, char **argv) {

    Eigen::VectorXd x = Eigen::VectorXd::Ones(1e5);
    x = 0.5*x;
    Eigen::VectorXd y = x;
    clock_t begin_time = clock();
    y = dnorm(x);
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    begin_time = clock();
    y = dnorm_gsl(x);
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    begin_time = clock();
    y = pnorm(x);
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    begin_time = clock();
    y = pnorm_gsl(x);
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    begin_time = clock();
    y = qnorm(x);
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    begin_time = clock();
    y = qnorm_gsl(x);
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

    return 0;
}
