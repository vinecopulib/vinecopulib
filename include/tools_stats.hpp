// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <random>
#include "tools_eigen.hpp"
#include "tools_c.h"

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
inline MatXd simulate_uniform(int n, int d)
{
    if ((n < 1) | (d < 1))
        throw std::runtime_error("both n and d must be at least 1.");
    MatXd U(n, d);
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return U.unaryExpr([&](double) { return distribution(generator); });
}

inline double pairwise_ktau(MatXd& u)
{
    double tau;
    int n = u.rows();
    int two = 2;
    ktau_matrix(u.data(), &two, &n, &tau);
    return tau;
}

inline double pairwise_cor(const MatXd& z)
{
    double rho;
    MatXd x = z.rowwise() - z.colwise().mean();
    MatXd sigma = x.adjoint() * x;
    rho = sigma(1,0) / sqrt(sigma(0,0) * sigma(1,1));

    return rho;
}
