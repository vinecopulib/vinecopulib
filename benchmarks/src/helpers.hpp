// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <string>
#include <vector>
#include <vinecopulib.hpp>

namespace bench {

//! Representative parameters per family (mid-range dependence, valid for
//! all rotations used below).
inline Eigen::MatrixXd
family_parameters(vinecopulib::BicopFamily family)
{
  using namespace vinecopulib;
  Eigen::MatrixXd par;
  switch (family) {
    case BicopFamily::indep:
      par.resize(0, 1);
      break;
    case BicopFamily::gaussian:
      par = Eigen::MatrixXd::Constant(1, 1, 0.5);
      break;
    case BicopFamily::student:
      par.resize(2, 1);
      par << 0.5, 4.0;
      break;
    case BicopFamily::clayton:
      par = Eigen::MatrixXd::Constant(1, 1, 3.0);
      break;
    case BicopFamily::gumbel:
      par = Eigen::MatrixXd::Constant(1, 1, 2.0);
      break;
    case BicopFamily::frank:
      par = Eigen::MatrixXd::Constant(1, 1, 5.0);
      break;
    case BicopFamily::joe:
      par = Eigen::MatrixXd::Constant(1, 1, 3.0);
      break;
    case BicopFamily::bb1:
      par.resize(2, 1);
      par << 1.0, 1.5;
      break;
    case BicopFamily::bb6:
      par.resize(2, 1);
      par << 1.5, 1.5;
      break;
    case BicopFamily::bb7:
      par.resize(2, 1);
      par << 1.5, 1.0;
      break;
    case BicopFamily::bb8:
      par.resize(2, 1);
      par << 3.0, 0.8;
      break;
    case BicopFamily::tawn:
      // (psi1, psi2, theta)
      par.resize(3, 1);
      par << 0.8, 0.5, 2.0;
      break;
    default:
      throw std::runtime_error("no benchmark parameters for this family");
  }
  return par;
}

//! Simulates copula data from a parametric model (deterministic via seed).
inline Eigen::MatrixXd
sim_data(vinecopulib::BicopFamily family, int rotation, size_t n, int seed = 5)
{
  vinecopulib::Bicop bc(family, rotation, family_parameters(family));
  return bc.simulate(n, false, { seed });
}

//! Discretizes the first column of a 2-column matrix to `levels` support
//! points and returns the 4-column (u, u^-) format used for var_types
//! {"d", "c"}.
inline Eigen::MatrixXd
discretize_first(const Eigen::MatrixXd& u, double levels = 10.0)
{
  Eigen::MatrixXd u_disc(u.rows(), 4);
  u_disc.col(0) = (u.col(0).array() * levels).ceil() / levels;
  u_disc.col(1) = u.col(1);
  u_disc.col(2) = (u.col(0).array() * levels).floor() / levels;
  u_disc.col(3) = u.col(1);
  return u_disc;
}

//! Simulates data from a Gaussian d-dimensional vine (deterministic).
inline vinecopulib::Vinecop
make_gaussian_vine(size_t d, int seed = 5)
{
  using namespace vinecopulib;
  auto structure = RVineStructure::simulate(d, false, { seed });
  auto pc_store = Vinecop::make_pair_copula_store(d);
  for (auto& tree : pc_store) {
    for (auto& pc : tree) {
      pc =
        Bicop(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.5));
    }
  }
  return Vinecop(structure, pc_store);
}

//! Simulates from a d-dimensional vine whose pair-copulas are all `family`
//! with parameter `par` (deterministic).
inline vinecopulib::Vinecop
make_vine(vinecopulib::BicopFamily family, double par, size_t d, int seed = 5)
{
  using namespace vinecopulib;
  auto structure = RVineStructure::simulate(d, false, { seed });
  auto pc_store = Vinecop::make_pair_copula_store(d);
  for (auto& tree : pc_store) {
    for (auto& pc : tree) {
      pc = Bicop(family, 0, Eigen::MatrixXd::Constant(1, 1, par));
    }
  }
  return Vinecop(structure, pc_store);
}

} // namespace bench
