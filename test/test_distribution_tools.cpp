/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "include/distribution_tools.hpp"
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
    time("Boost pnom:", [m]{ auto a = pnorm(m); });
    time("GSL pnorm:  ", [m]{ auto b = pnorm_gsl(m); });
    time("Boost qnom:", [m]{ auto a = qnorm(m); });
    time("GSL qnorm:  ", [m]{ auto b = qnorm_gsl(m); });

    return 0;
}
