// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
namespace tools_select {
//! @brief Creates the candidate models for a pair: the eligible families
//! with the rotations that fit the association direction.
//! @param data Captured by reference to avoid data copies;
//!     should NOT be modified though.
//! @param controls The controls for the pair-copula fit, whose family set and
//!     selection criterion the candidates are drawn from.
//! @param var_types The pair's variable types; they decide which families
//!     are eligible and how rotations are chosen.
inline std::vector<Bicop>
create_candidate_bicops(const Eigen::MatrixXd& data,
                        const FitControlsBicop& controls,
                        const std::vector<std::string>& var_types)
{
  std::vector<BicopFamily> families =
    get_candidate_families(controls, var_types);

  // check whether dependence is negative or positive
  double tau = wdm::wdm(data.leftCols(2), "tau", controls.get_weights())(0, 1);
  std::vector<int> which_rotations;
  if (controls.get_allow_rotations()) {
    if (tau > 0) {
      which_rotations = { 0, 180 };
    } else {
      which_rotations = { 90, 270 };
    }
  }

  // create Bicop objects for all valid family/rotation combinations, each
  // carrying the pair's variable types
  std::vector<Bicop> new_bicops;
  const Eigen::MatrixXd no_pars;
  auto add = [&](BicopFamily fam, int rotation) {
    new_bicops.emplace_back(fam, rotation, no_pars, var_types);
  };
  for (auto& fam : families) {
    if (tools_stl::is_member(fam, bicop_families::rotationless)) {
      add(fam, 0);
    } else if (tools_stl::is_member(fam, bicop_families::two_rotations)) {
      // both circular orientations; linear tau carries no information here
      add(fam, 0);
      if (controls.get_allow_rotations()) {
        add(fam, 90);
      }
    } else {
      if (controls.get_allow_rotations()) {
        add(fam, which_rotations[0]);
        add(fam, which_rotations[1]);
      } else if (tau > 0) {
        add(fam, 0);
      }
    }
  }

  // remove combinations based on symmetry characteristics
  if (controls.get_preselect_families()) {
    preselect_candidates(
      new_bicops, data.leftCols(2), tau, controls.get_weights());
  }

  return new_bicops;
}

//! @brief The families to try for a pair: the family set of the controls
//! (all families when empty), restricted to the parametric method and to the
//! families eligible for the pair's variable types.
//! @param controls The fit controls.
//! @param var_types The pair's variable types.
inline std::vector<BicopFamily>
get_candidate_families(const FitControlsBicop& controls,
                       const std::vector<std::string>& var_types)
{
  std::vector<BicopFamily> family_set = controls.get_family_set();
  if (family_set.empty()) {
    // use all (allowed) families
    if (controls.get_parametric_method() == "itau") {
      family_set = bicop_families::itau;
    } else {
      family_set = bicop_families::all;
    }
  } else {
    if (controls.get_parametric_method() == "itau") {
      family_set = tools_stl::intersect(family_set, bicop_families::itau);
      if (family_set.empty()) {
        throw std::runtime_error("No family with method itau provided");
      }
    }
  }

  family_set = eligible_families(family_set, var_types);
  if (family_set.empty()) {
    std::string eligible;
    for (auto fam : eligible_families(bicop_families::all, var_types)) {
      eligible += (eligible.empty() ? "" : ", ") + get_family_name(fam);
    }
    throw std::runtime_error(
      "no family in the family set can model a pair with var_types (" +
      var_types[0] + ", " + var_types[1] +
      "); eligible families are: " + eligible);
  }

  return family_set;
}

//! removes candidates whose symmetry properties does not correspond to those
//! of the data.
inline void
preselect_candidates(std::vector<Bicop>& bicops,
                     const Eigen::MatrixXd& data,
                     double tau,
                     const Eigen::VectorXd& weights)
{
  auto c = get_c1c2(data, tau, weights);
  bicops.erase(std::remove_if(bicops.begin(),
                              bicops.end(),
                              [&](const Bicop& cop) {
                                return !(preselect_family(c, tau, cop));
                              }),
               bicops.end());
}

inline std::vector<double>
get_c1c2(const Eigen::MatrixXd& data,
         double tau,
         const Eigen::VectorXd& weights)
{
  size_t n = data.rows();
  Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, 2);
  Eigen::MatrixXd z1 = x;
  Eigen::MatrixXd z2 = x;
  x = tools_stats::qnorm(data);

  int count1 = 0, count2 = 0;
  for (size_t j = 0; j < n; ++j) {
    if (tau > 0) {
      if ((x(j, 0) > 0) && (x(j, 1) > 0)) {
        z1.row(count1) = x.row(j);
        ++count1;
      }
      if ((x(j, 0) < 0) && (x(j, 1) < 0)) {
        z2.row(count2) = x.row(j);
        ++count2;
      }
    } else {
      if ((x(j, 0) < 0) && (x(j, 1) > 0)) {
        z1.row(count1) = x.row(j);
        ++count1;
      }
      if ((x(j, 0) > 0) && (x(j, 1) < 0)) {
        z2.row(count2) = x.row(j);
        ++count2;
      }
    }
  }

  // if one of the quadrants is empty, we see it as independent
  double c1, c2;
  Eigen::VectorXd w;

  if (count1 == 0) {
    c1 = 0.0;
  } else {
    w = (weights.size() > 0) ? weights.head(count1 - 1) : weights;
    c1 = wdm::wdm(z1.block(0, 0, count1 - 1, 2), "cor", w)(0, 1);
  }

  if (count2 == 0) {
    c2 = 0.0;
  } else {
    w = (weights.size() > 0) ? weights.head(count2 - 1) : weights;
    c2 = wdm::wdm(z2.block(0, 0, count2 - 1, 2), "cor", w)(0, 1);
  }

  return { c1, c2 };
}

inline bool
preselect_family(std::vector<double> c, double tau, const Bicop& bicop)
{
  using namespace tools_stl;
  BicopFamily family = bicop.get_family();
  int rotation = bicop.get_rotation();

  bool preselect = false;
  if (is_member(family, bicop_families::circular)) {
    // the tail and sign heuristics below are defined for linear dependence
    preselect = true;
  } else if (is_member(family, bicop_families::rotationless)) {
    preselect = true;
    if ((std::fabs(c[0] - c[1]) > 0.3) & (family == BicopFamily::frank))
      preselect = false;
  } else {
    if (is_member(family, bicop_families::bb)) {
      if ((tau > 0) && is_member(rotation, { 0, 180 })) {
        preselect = true;
      }
      if ((tau < 0) && is_member(rotation, { 90, 270 })) {
        preselect = true;
      }
    }
    bool is_90or180 = is_member(rotation, { 90, 180 });
    if (c[0] - c[1] > 0.05) {
      if (is_member(family, bicop_families::lt) && is_90or180) {
        preselect = true;
      }
      if (is_member(family, bicop_families::ut) && !is_90or180) {
        preselect = true;
      }
    } else if (c[0] - c[1] < -0.05) {
      if (is_member(family, bicop_families::lt) && !is_90or180) {
        preselect = true;
      }
      if (is_member(family, bicop_families::ut) && is_90or180) {
        preselect = true;
      }
    } else {
      if ((tau > 0) && is_member(rotation, { 0, 180 })) {
        preselect = true;
      }
      if ((tau < 0) && is_member(rotation, { 90, 270 })) {
        preselect = true;
      }
    }
  }
  return preselect;
}
}
}
