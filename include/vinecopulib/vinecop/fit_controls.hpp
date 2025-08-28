// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vinecopulib/bicop/fit_controls.hpp>

namespace vinecopulib {
//! @brief A class for controlling fits of vine copula models.
//!
class FitControlsVinecop : public FitControlsBicop
{
public:
  // Constructor
  [[deprecated("Use FitControls instead")]] FitControlsVinecop();

  [[deprecated("Use FitControls instead")]] explicit FitControlsVinecop(
    std::vector<BicopFamily> family_set,
    std::string parametric_method = "mle",
    std::string nonparametric_method = "constant",
    double nonparametric_mult = 1.0,
    size_t nonparametric_grid_size = 30,
    size_t trunc_lvl = std::numeric_limits<size_t>::max(),
    std::string tree_criterion = "tau",
    double threshold = 0.0,
    std::string selection_criterion = "aic",
    const Eigen::VectorXd& weights = Eigen::VectorXd(),
    double psi0 = 0.9,
    bool preselect_families = true,
    bool select_trunc_lvl = false,
    bool select_threshold = false,
    bool select_families = true,
    bool show_trace = false,
    size_t num_threads = 1,
    std::string tree_algorithm = "mst_prim",
    bool allow_rotations = true,
    std::vector<int> seeds = std::vector<int>());

  [[deprecated("Use FitControls instead")]] explicit FitControlsVinecop(
    const FitControlsBicop& controls,
    size_t trunc_lvl = std::numeric_limits<size_t>::max(),
    std::string tree_criterion = "tau",
    double threshold = 0.0,
    bool select_trunc_lvl = false,
    bool select_threshold = false,
    bool select_families = true,
    bool show_trace = false,
    std::string tree_algorithm = "mst_prim",
    std::vector<int> seeds = std::vector<int>());

  explicit FitControlsVinecop(const FitControls& controls);

  [[deprecated("Use get_trunc_lvl() instead")]] size_t get_truncation_level()
    const;
  size_t get_trunc_lvl() const;

  [[deprecated("Use set_trunc_lvl() instead")]] void set_truncation_level(
    size_t trunc_lvl);
  void set_trunc_lvl(size_t trunc_lvl);

  std::string get_tree_criterion() const;
  void set_tree_criterion(std::string tree_criterion);

  double get_threshold() const;
  void set_threshold(double threshold);

  bool get_show_trace() const;
  void set_show_trace(bool show_trace);

  [[deprecated("Use get_select_trunc_lvl() instead")]] bool
  get_select_truncation_level() const;
  bool get_select_trunc_lvl() const;
  [[deprecated("Use set_select_trunc_lvl() instead")]] void
  set_select_truncation_level(bool select_trunc_lvl);
  void set_select_trunc_lvl(bool select_trunc_lvl);

  bool get_select_threshold() const;
  void set_select_threshold(bool select_threshold);

  bool get_select_families() const;
  void set_select_families(bool select_families);

  std::string get_tree_algorithm() const;
  void set_tree_algorithm(std::string tree_algorithm);

  std::vector<int> get_seeds() const;
  void set_seeds(std::vector<int> seeds);

  FitControlsBicop get_fit_controls_bicop() const;
  void set_fit_controls_bicop(FitControlsBicop controls);

  // Misc
  bool needs_sparse_select() const;
  std::string str() const;
};
}

#include <vinecopulib/vinecop/implementation/fit_controls.ipp>
