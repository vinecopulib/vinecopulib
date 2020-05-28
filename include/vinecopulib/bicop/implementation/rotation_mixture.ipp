// Copyright © 2016-2020 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/clayton.hpp>
#include <vinecopulib/bicop/gumbel.hpp>
#include <vinecopulib/bicop/joe.hpp>
#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline RotationMixtureBicop::RotationMixtureBicop()
{
  base_bicop_ = BicopPtr(new IndepBicop());
  parameters_ = Eigen::VectorXd::Zero(4);
  parameters_lower_bounds_ = Eigen::VectorXd::Zero(4);
  parameters_upper_bounds_ = Eigen::VectorXd::Ones(4);
}

inline Eigen::VectorXd
RotationMixtureBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd pdf = Eigen::VectorXd::Zero(u.rows());
  for (int i = 0; i < 4; i++) {
    base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(i)));
    pdf += base_bicop_->pdf(uu) / 4;
    rotate_90(uu);
  }
  return pdf;
}

inline Eigen::VectorXd
RotationMixtureBicop::cdf(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd cdf = Eigen::VectorXd::Zero(u.rows());
  for (int i = 0; i < 4; i++) {
    base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(i)));
    cdf += base_bicop_->cdf(uu) / 4;
    rotate_90(uu);
  }
  return cdf;
}

inline Eigen::VectorXd
RotationMixtureBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd h1 = Eigen::VectorXd::Zero(u.rows());
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(u.rows());

  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(0)));
  h1 += 0.25 * base_bicop_->hfunc1(uu);

  rotate_90(uu);
  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(1)));
  h1 += 0.25 * base_bicop_->hfunc2(uu);

  rotate_90(uu);
  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(2)));
  h1 += 0.25 * (ones - base_bicop_->hfunc1(uu));

  rotate_90(uu);
  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(3)));
  h1 += 0.25 * (ones - base_bicop_->hfunc2(uu));

  return h1;
}

inline Eigen::VectorXd
RotationMixtureBicop::hfunc2_raw(const Eigen::MatrixXd& u)
{
  auto uu = u;
  Eigen::VectorXd h2 = Eigen::VectorXd::Zero(u.rows());
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(u.rows());

  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(0)));
  h2 += 0.25 * base_bicop_->hfunc2(uu);

  rotate_90(uu);
  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(1)));
  h2 += 0.25 * (ones - base_bicop_->hfunc1(uu));

  rotate_90(uu);
  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(2)));
  h2 += 0.25 * (ones - base_bicop_->hfunc2(uu));

  rotate_90(uu);
  base_bicop_->set_parameters(Eigen::VectorXd::Constant(1, parameters_(3)));
  h2 += 0.25 * base_bicop_->hfunc1(uu);

  return h2;
}

inline Eigen::VectorXd
RotationMixtureBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  return hinv1_num(u);
}

inline Eigen::VectorXd
RotationMixtureBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  return hinv2_num(u);
}

inline Eigen::MatrixXd 
RotationMixtureBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

inline double
RotationMixtureBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  auto oldpars = this->get_parameters();
  auto old_types = var_types_;
  this->set_parameters(parameters);
  var_types_ = { "c", "c" };

  std::vector<int> seeds = {
    204967043, 733593603, 184618802, 399707801, 290266245
  };
  auto u = tools_stats::ghalton(1000, 2, seeds);
  u.col(1) = hinv1_raw(u);

  this->set_parameters(oldpars);
  var_types_ = old_types;
  return wdm::wdm(u, "tau")(0, 1);
}

//! @brief Rotates the data corresponding to the models rotation.
//! @param u An `n x 2` matrix.
inline void
RotationMixtureBicop::rotate_90(Eigen::MatrixXd& u) const
{
  // counter-clockwise rotations
  u.col(0).swap(u.col(1));
  u.col(1) = 1 - u.col(1).array();
  if (u.cols() == 4) {
    u.col(2).swap(u.col(3));
    u.col(3) = 1 - u.col(3).array();
  }
}

// Gumbel mixture --------------------

inline RMGumbelBicop::RMGumbelBicop()
{
  family_ = BicopFamily::rmgumbel;
  parameters_ = Eigen::VectorXd::Constant(4, 1.0);
  parameters_lower_bounds_ = Eigen::VectorXd::Ones(4);
  parameters_upper_bounds_ = Eigen::VectorXd::Ones(4) * 50;
  base_bicop_ = BicopPtr(new GumbelBicop());
}

inline Eigen::VectorXd
RMGumbelBicop::get_start_parameters(const double)
{
  return Eigen::VectorXd::Ones(4) * 2.0;
}


// Clayton mixture --------------------

inline RMClaytonBicop::RMClaytonBicop()
{
  family_ = BicopFamily::rmclayton;
  parameters_ = Eigen::VectorXd::Constant(4, 1e-10);
  parameters_lower_bounds_ = Eigen::VectorXd::Ones(4) * 1e-10;
  parameters_upper_bounds_ = Eigen::VectorXd::Ones(4) * 28;
  base_bicop_ = BicopPtr(new ClaytonBicop());
}

inline Eigen::VectorXd
RMClaytonBicop::get_start_parameters(const double)
{
  return Eigen::VectorXd::Ones(4) * 1.5;
}


// Joe mixture --------------------

inline RMJoeBicop::RMJoeBicop()
{
  family_ = BicopFamily::rmjoe;
  parameters_ = Eigen::VectorXd::Constant(4, 1.0);
  parameters_lower_bounds_ = Eigen::VectorXd::Ones(4);
  parameters_upper_bounds_ = Eigen::VectorXd::Ones(4) * 30;
  base_bicop_ = BicopPtr(new JoeBicop());
}

inline Eigen::VectorXd
RMJoeBicop::get_start_parameters(const double)
{
  return Eigen::VectorXd::Ones(4) * 1.5;
}

}
