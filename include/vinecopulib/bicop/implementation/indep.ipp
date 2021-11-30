// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

inline IndepBicop::IndepBicop()
{
  family_ = BicopFamily::indep;
  parameters_ = Matrix();
}

inline Vector
IndepBicop::pdf_raw(const Matrix& u)
{
  auto f = [](double, double) { return 1.0; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Vector
IndepBicop::cdf(const Matrix& u)
{
  return u.rowwise().prod();
}

inline Vector
IndepBicop::hfunc1_raw(const Matrix& u)
{
  auto f = [](double, double u2) { return u2; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Vector
IndepBicop::hfunc2_raw(const Matrix& u)
{
  auto f = [](double u1, double) { return u1; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Vector
IndepBicop::hinv1_raw(const Matrix& u)
{
  auto f = [](double, double u2) { return u2; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Vector
IndepBicop::hinv2_raw(const Matrix& u)
{
  auto f = [](double u1, double) { return u1; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Matrix
IndepBicop::tau_to_parameters(const double&)
{
  return Vector();
}

inline double
IndepBicop::parameters_to_tau(const Matrix&)
{
  return 0.0;
}

inline Vector
IndepBicop::get_start_parameters(const double tau)
{
  return tau_to_parameters(tau);
}

inline void
IndepBicop::flip()
{
  // nothing to do because independence copula is radially syemmetric
}
}
