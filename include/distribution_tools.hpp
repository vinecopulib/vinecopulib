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

#ifndef VINECOPULIB_DISTRIBUTION_TOOLS_HPP
#define VINECOPULIB_DISTRIBUTION_TOOLS_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <Eigen/Dense>
#include <random>

template<typename T> T dnorm(const T& x)
{
    boost::math::normal std_normal;
    return x.unaryExpr(boost::bind<double>(boost::math::pdf<boost::math::normal,double>, std_normal, _1));
};

template<typename T> T pnorm(const T& x)
{
    boost::math::normal std_normal;
    return x.unaryExpr(boost::bind<double>(boost::math::cdf<boost::math::normal,double>, std_normal, _1));
};

template<typename T> T qnorm(const T& x)
{
    boost::math::normal std_normal;
    return x.unaryExpr(boost::bind<double>(boost::math::quantile<boost::math::normal,double>, std_normal, _1));
};

template<typename T> T dt(const T& x, double nu)
{
    boost::math::students_t dist(nu);
    return x.unaryExpr(boost::bind<double>(boost::math::pdf<boost::math::students_t,double>, dist, _1));
};

template<typename T> T pt(const T& x, double nu)
{
    boost::math::students_t dist(nu);
    return x.unaryExpr(boost::bind<double>(boost::math::cdf<boost::math::students_t,double>, dist, _1));
};

template<typename T> T qt(const T& x, double nu)
{
    boost::math::students_t dist(nu);
    return x.unaryExpr(boost::bind<double>(boost::math::quantile<boost::math::students_t,double>, dist, _1));
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

#endif //VINECOPULIB_DISTRIBUTION_TOOLS_HPP
