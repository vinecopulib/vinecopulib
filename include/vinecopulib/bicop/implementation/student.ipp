// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
inline StudentBicop::StudentBicop()
{
  family_ = BicopFamily::student;
  parameters_ = Vector(2);
  parameters_lower_bounds_ = Vector(2);
  parameters_upper_bounds_ = Vector(2);
  parameters_ << 0, 50;
  parameters_lower_bounds_ << -1, 2;
  parameters_upper_bounds_ << 1, 50;
}

inline Vector
StudentBicop::pdf_raw(const Matrix& u)
{
  double rho = double(this->parameters_(0));
  double nu = double(this->parameters_(1));
  Vector f = Vector::Ones(u.rows());
  Matrix tmp = tools_stats::qt(u, nu);

  f = tmp.col(0).cwiseAbs2() + tmp.col(1).cwiseAbs2() -
      (2 * rho) * tmp.rowwise().prod();
  f /= nu * (1.0 - pow(rho, 2.0));
  f = f + Vector::Ones(u.rows());
  f = f.array().pow(-(nu + 2.0) / 2.0);
  f = f.cwiseQuotient(tools_stats::dt(tmp, nu).rowwise().prod());
  f *= boost::math::tgamma_ratio((nu + 2.0) / 2.0, nu / 2.0);
  f /= (nu * constant::pi * sqrt(1.0 - pow(rho, 2.0)));

  return f;
}

inline Vector
StudentBicop::cdf(const Matrix& u)
{
  using namespace tools_stats;

  double rho = double(this->parameters_(0));
  double nu = double(this->parameters_(1));

  // for integer nu, just use pbvt
  // otherwise, interpolate linearly between floor(nu) and ceil(nu)
  if (nu == round(nu)) {
    int inu = static_cast<int>(nu);
    return pbvt(qt(u, inu), inu, rho);
  } else {
    int nu1 = static_cast<int>(std::floor(nu));
    int nu2 = static_cast<int>(std::ceil(nu));
    double weight = (nu - static_cast<double>(nu1)) /
                    (static_cast<double>(nu2) - static_cast<double>(nu1));
    return pbvt(qt(u, nu1), nu1, rho) * (1 - weight) +
           pbvt(qt(u, nu2), nu2, rho) * weight;
  }
}

inline Vector
StudentBicop::hfunc1_raw(const Matrix& u)
{
  double rho = double(this->parameters_(0));
  double nu = double(this->parameters_(1));
  Vector h = Vector::Ones(u.rows());
  Matrix tmp = tools_stats::qt(u, nu);
  h = nu * h + tmp.col(0).cwiseAbs2();
  h *= (1.0 - pow(rho, 2)) / (nu + 1.0);
  h = h.cwiseSqrt().cwiseInverse().cwiseProduct(tmp.col(1) - rho * tmp.col(0));
  h = tools_stats::pt(h, nu + 1.0);

  return h;
}

inline Vector
StudentBicop::hinv1_raw(const Matrix& u)
{
  double rho = double(this->parameters_(0));
  double nu = double(this->parameters_(1));
  Vector hinv = Vector::Ones(u.rows());
  Vector tmp = u.col(1);
  Vector tmp2 = u.col(0);
  tmp = tools_stats::qt(tmp, nu + 1.0);
  tmp2 = tools_stats::qt(tmp2, nu);

  hinv = nu * hinv + tmp2.cwiseAbs2();
  hinv *= (1.0 - pow(rho, 2)) / (nu + 1.0);
  hinv = hinv.cwiseSqrt().cwiseProduct(tmp) + rho * tmp2;
  hinv = tools_stats::pt(hinv, nu);

  return hinv;
}

inline Vector
StudentBicop::get_start_parameters(const double tau)
{
  Vector parameters = get_parameters();
  parameters(0) = std::sin(tau * constant::pi / 2);
  parameters(1) = 5;
  return parameters;
}

inline Matrix
StudentBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}
}
