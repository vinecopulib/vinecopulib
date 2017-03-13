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

#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions.hpp>
#include <boost/function.hpp>
#include <Eigen/Dense>
#include <random>

template<typename T> T dnorm(const T& x)
{
    boost::math::normal dist;
    return x.unaryExpr([&dist](double y) {return boost::math::pdf(dist, y);});
};

template<typename T> T pnorm(const T& x)
{
    boost::math::normal dist;
    return x.unaryExpr([&dist](double y) {return boost::math::cdf(dist, y);});
};

template<typename T> T qnorm(const T& x)
{
    boost::math::normal dist;
    return x.unaryExpr([&dist](double y) {return boost::math::quantile(dist, y);});
};

template<typename T> T dt(const T& x, double nu)
{
    boost::math::students_t dist(nu);
    return x.unaryExpr([&dist](double y) {return boost::math::pdf(dist, y);});
};

template<typename T> T pt(const T& x, double nu)
{
    boost::math::students_t dist(nu);
    return x.unaryExpr([&dist](double y) {return boost::math::cdf(dist, y);});
};

template<typename T> T qt(const T& x, double nu)
{
    boost::math::students_t dist(nu);
    return x.unaryExpr([&dist](double y) {return boost::math::quantile(dist, y);});
};

/* A GSL ALTERNATIVE
 *
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <functional>

template<typename T> T dnorm(const T& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_ran_ugaussian_pdf));
};
template<typename T> T pnorm(const T& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_cdf_ugaussian_P));
};
template<typename T> T qnorm(const T& x)
{
    return x.unaryExpr(std::ptr_fun(gsl_cdf_ugaussian_Pinv));
};

template<typename T> T dt(const T& x, double nu)
{
    auto tpdf = std::bind(gsl_ran_tdist_pdf, std::placeholders::_1, nu);
    return x.unaryExpr(std::function<double(double)>(tpdf));
};
template<typename T> T pt(const T& x, double nu)
{
    auto tcdf = std::bind(gsl_cdf_tdist_P, std::placeholders::_1, nu);
    return x.unaryExpr(std::function<double(double)>(tcdf));
};
template<typename T> T qt(const T& x, double nu)
{
    auto tquantile = std::bind(gsl_cdf_tdist_Pinv, std::placeholders::_1, nu);
    return x.unaryExpr(std::function<double(double)>(tquantile));
};*/

//! Simulate from the multivariate uniform distribution
//! 
//! @param n number of observations.
//! @param d dimension.
//! 
//! @return A nxd matrix of independent U[0, 1] random variables.
inline Eigen::MatrixXd simulate_uniform(int n, int d) 
{
    if ((n < 1) | (d < 1))
        throw std::runtime_error("both n and d must be at least 1.");
    Eigen::MatrixXd U(n, d);
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return U.unaryExpr([&](double) { return distribution(generator); });
}

