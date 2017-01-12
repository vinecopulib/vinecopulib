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

#ifndef VINECOPULIB_BOOST_TOOLS_HPP
#define VINECOPULIB_BOOST_TOOLS_HPP

#include <Eigen/Dense>
typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

#include <boost/bind.hpp>
#include <boost/math/distributions.hpp>
#include <boost/function.hpp>

namespace bm = boost::math;

namespace distributions{
    bm::normal std_normal;
}

boost::function<double(double)> dnorm_scalar = boost::bind<double>(bm::pdf<bm::normal,double>,
                                                                   distributions::std_normal, _1);

boost::function<double(double)> pnorm_scalar = boost::bind<double>(bm::cdf<bm::normal,double>,
                                                                   distributions::std_normal, _1);

boost::function<double(double)> qnorm_scalar = boost::bind<double>(bm::quantile<bm::normal,double>,
                                                                   distributions::std_normal, _1);


VecXd dnorm(const VecXd& x) { return x.unaryExpr(dnorm_scalar); }

VecXd pnorm(const VecXd& x) { return x.unaryExpr(pnorm_scalar); }

VecXd qnorm(const VecXd& p) { return p.unaryExpr(qnorm_scalar); }

//namespace distributions{
//    bm::normal std_normal(0,1);
//}

//bm::normal std_normal(0,1);
//boost::function<double(double)> qnorm = boost::bind(quantile, distributions::std_normal, _1);
//boost::function<double(double)> qnorm = boost::bind(bm::cdf<double>, boost::ref(std_normal), _1);
//VecXd qsnorm(VecXd u);

#endif //VINECOPULIB_BOOST_TOOLS_HPP
