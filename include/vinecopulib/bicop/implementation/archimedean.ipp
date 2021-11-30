// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
// inline Vector ArchimedeanBicop::pdf(
//    const Matrix &u
//)
//{
//    auto f = [this](const double &u1, const double &u2) {
//        double temp = generator_inv(generator(u1) + generator(u2));
//        temp = log(std::abs(generator_derivative2(temp))) -
//            3.0 * log(std::abs(generator_derivative(temp)));
//        temp += std::log(std::abs(generator_derivative(u1)));
//        temp += std::log(std::abs(generator_derivative(u2)));
//        return std::exp(temp);
//    };
//    return tools_eigen::binaryExpr_or_nan(u, f);
//}

inline Vector
ArchimedeanBicop::cdf(const Matrix& u)
{
  auto f = [this](const double& u1, const double& u2) {
    return generator_inv(generator(u1) + generator(u2));
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Vector
ArchimedeanBicop::hfunc1_raw(const Matrix& u)
{
  auto f = [this](const double& u1, const double& u2) {
    double temp = generator_inv(generator(u1) + generator(u2));
    temp = generator_derivative(u1) / generator_derivative(temp);
    return std::isnan(temp) ? u2 : std::min(temp, 1.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Vector
ArchimedeanBicop::hfunc2_raw(const Matrix& u)
{
  return hfunc1_raw(tools_eigen::swap_cols(u));
}

inline Vector
ArchimedeanBicop::hinv1_raw(const Matrix& u)
{
  Vector hinv = hinv1_num(u);
  return hinv;
}

inline Vector
ArchimedeanBicop::hinv2_raw(const Matrix& u)
{
  return hinv1_raw(tools_eigen::swap_cols(u));
}

inline Vector
ArchimedeanBicop::get_start_parameters(const double)
{
  Matrix lb = this->get_parameters_lower_bounds();
  Vector parameters = lb + Vector::Constant(2, 0.1);
  return parameters;
}
}
