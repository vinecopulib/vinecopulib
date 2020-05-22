// Copyright © 2016-2020 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/rmgumbel.hpp>

namespace vinecopulib {
inline RMGumbelBicop::RMGumbelBicop()
{
  family_ = BicopFamily::rmgumbel;
  parameters_ = Eigen::VectorXd::Constant(4, 1.0);
  parameters_lower_bounds_ = Eigen::VectorXd::Ones(4);
  parameters_upper_bounds_ = Eigen::VectorXd::Ones(4) * 50;
}

inline Eigen::VectorXd
RMGumbelBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd pdf = Eigen::VectorXd::Zero(u.rows());
  GumbelBicop bc;
  for (int i = 0; i < 4; i++) {
    bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(i)));
    pdf += bc.pdf(uu) / 4;
    rotate_90(uu);
  }
  return pdf;
}

inline Eigen::VectorXd
RMGumbelBicop::cdf(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd cdf = Eigen::VectorXd::Zero(u.rows());
  GumbelBicop bc;
  for (int i = 0; i < 4; i++) {
    bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(i)));
    cdf += bc.cdf(uu) / 4;
    rotate_90(uu);
  }
  return cdf;
}

inline Eigen::VectorXd
RMGumbelBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd h1 = Eigen::VectorXd::Zero(u.rows());
  GumbelBicop bc;
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(u.rows());

  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(0)));
  h1 += 0.25 * bc.hfunc1(uu);

  rotate_90(uu);
  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(1)));
  h1 += 0.25 * bc.hfunc2(uu);

  rotate_90(uu);
  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(2)));
  h1 += 0.25 * (ones - bc.hfunc1(uu));

  rotate_90(uu);
  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(3)));
  h1 += 0.25 * (ones - bc.hfunc2(uu));

  return h1;
}

inline Eigen::VectorXd
RMGumbelBicop::hfunc2_raw(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd h2 = Eigen::VectorXd::Zero(u.rows());
  GumbelBicop bc;
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(u.rows());

  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(0)));
  h2 += 0.25 * bc.hfunc2(uu);

  rotate_90(uu);
  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(1)));
  h2 += 0.25 * (ones - bc.hfunc1(uu));

  rotate_90(uu);
  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(2)));
  h2 += 0.25 * (ones - bc.hfunc2(uu));

  rotate_90(uu);
  bc.set_parameters(Eigen::VectorXd::Constant(1, parameters_(3)));
  h2 += 0.25 * bc.hfunc1(uu);

  return h2;
}

inline Eigen::VectorXd
RMGumbelBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  return hinv1_num(u);
}

inline Eigen::VectorXd
RMGumbelBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  return hinv2_num(u);
}

inline Eigen::VectorXd
RMGumbelBicop::get_start_parameters(const double tau)
{
  return Eigen::VectorXd::Ones(4) * 2.0;
}

inline Eigen::MatrixXd 
RMGumbelBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

inline double
RMGumbelBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  return 0.0;
}


//! @brief Rotates the data corresponding to the models rotation.
//! @param u An `n x 2` matrix.
inline void
RMGumbelBicop::rotate_90(Eigen::MatrixXd& u) const
{
  // counter-clockwise rotations
  u.col(0).swap(u.col(1));
  u.col(1) = 1 - u.col(1).array();
  if (u.cols() == 4) {
    u.col(2).swap(u.col(3));
    u.col(3) = 1 - u.col(3).array();
  }
}
}
