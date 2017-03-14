// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "tools_eigen.hpp"

//! Numerical inversion of a function
//!
//! Computes the inverse \f$f^{-1}\f$ of a function \f$f\f$ by the bisection
//! method.
//!
//! @param x evaluation points.
//! @param f the function to invert.
//! @param lb lower bound.
//! @param ub upper bound.
//! @param n_iter the number of iterations for the bisection (defaults to 35,
//! guaranteeing an accuracy of 0.5^35 ~= 6e-11).
//!
//! @return f^{-1}(x).
VecXd invert_f(const VecXd& x, std::function<VecXd(const VecXd&)> f, const double lb, const double ub, int n_iter)
{
    VecXd xl = VecXd::Constant(x.size(), lb);
    VecXd xh = VecXd::Constant(x.size(), ub);
    VecXd x_tmp = x;
    for (int iter = 0; iter < n_iter; ++iter) {
        x_tmp = (xh + xl) / 2.0;
        VecXd fm = f(x_tmp) - x;
        xl = (fm.array() < 0).select(x_tmp, xl);
        xh = (fm.array() < 0).select(xh, x_tmp);
    }

    return x_tmp;
}
