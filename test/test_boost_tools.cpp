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

#include "include/boost_tools.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <Eigen/Dense>
#include <chrono>
#include <iostream>

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

template <typename T>
void time(const std::string label, const T &it)
{
    auto start = std::chrono::high_resolution_clock::now();
    it();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << label << "Elapsed time: " << elapsed.count() << " s\n";
}


int main(int argc, char **argv) {
    Eigen::MatrixXd m = Eigen::MatrixXd::Ones(10000, 10000);
    m *= 0.5;

    time("Boost dnorm:", [m]{ auto a = dnorm(m); });
    time("GSL dnorm:  ", [m]{ auto b = dnorm_gsl(m); });
    time("Boost pnorm:", [m]{ auto a = pnorm(m); });
    time("GSL pnorm:  ", [m]{ auto b = pnorm_gsl(m); });
    time("Boost qnom:", [m]{ auto a = qnorm(m); });
    time("GSL qnorm:  ", [m]{ auto b = qnorm_gsl(m); });

    return 0;
}
