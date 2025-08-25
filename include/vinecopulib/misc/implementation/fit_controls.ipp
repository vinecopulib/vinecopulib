// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <limits>
#include <stdexcept>
#include <thread>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

//! @brief Creates default configuration for bivariate copula models.
inline FitControlsConfig
FitControlsConfig::bicop_defaults()
{
  FitControlsConfig config;
  config.family_set = bicop_families::all;
  config.parametric_method = "mle";
  config.nonparametric_method = "constant";
  config.nonparametric_mult = 1.0;
  config.nonparametric_grid_size = 30;
  config.selection_criterion = "aic";
  config.weights = Eigen::VectorXd();
  config.psi0 = 0.9;
  config.preselect_families = true;
  config.allow_rotations = true;
  config.num_threads = 1;
  return config;
}

//! @brief Creates default configuration for vine copula models.
inline FitControlsConfig
FitControlsConfig::vinecop_defaults()
{
  FitControlsConfig config = bicop_defaults(); // Start with bicop defaults

  // Add vine-specific defaults
  config.trunc_lvl = std::numeric_limits<size_t>::max();
  config.tree_criterion = "tau";
  config.threshold = 0.0;
  config.select_trunc_lvl = false;
  config.select_threshold = false;
  config.select_families = true;
  config.show_trace = false;
  config.tree_algorithm = "mst_prim";
  config.seeds = std::vector<int>();

  return config;
}

//! @brief Validates all set values and applies defaults for unset values.
inline void
FitControlsConfig::validate_and_set_defaults()
{
  // Apply defaults for unset values, then validate
  if (!optional::has_value(family_set)) {
    family_set = bicop_families::all;
  }

  if (!optional::has_value(parametric_method)) {
    parametric_method = "mle";
  }
  check_parametric_method();

  if (!optional::has_value(nonparametric_method)) {
    nonparametric_method = "constant";
  }
  check_nonparametric_method();

  if (!optional::has_value(nonparametric_mult)) {
    nonparametric_mult = 1.0;
  }
  check_nonparametric_mult();

  if (!optional::has_value(nonparametric_grid_size)) {
    nonparametric_grid_size = 30;
  }
  check_nonparametric_grid_size();

  if (!optional::has_value(selection_criterion)) {
    selection_criterion = "aic";
  }
  check_selection_criterion();

  if (!optional::has_value(weights)) {
    weights = Eigen::VectorXd();
  }

  if (!optional::has_value(psi0)) {
    psi0 = 0.9;
  }
  check_psi0();

  if (!optional::has_value(preselect_families)) {
    preselect_families = true;
  }

  if (!optional::has_value(allow_rotations)) {
    allow_rotations = true;
  }

  if (!optional::has_value(num_threads)) {
    num_threads = 1;
  }
  // Process num_threads (validation and adjustment)
  num_threads = process_num_threads(optional::value(num_threads));

  // Vine-specific defaults and validation
  if (!optional::has_value(trunc_lvl)) {
    trunc_lvl = std::numeric_limits<size_t>::max();
  }

  if (!optional::has_value(tree_criterion)) {
    tree_criterion = "tau";
  }
  check_tree_criterion();

  if (!optional::has_value(threshold)) {
    threshold = 0.0;
  }
  check_threshold();

  if (!optional::has_value(select_trunc_lvl)) {
    select_trunc_lvl = false;
  }

  if (!optional::has_value(select_threshold)) {
    select_threshold = false;
  }

  if (!optional::has_value(select_families)) {
    select_families = true;
  }

  if (!optional::has_value(show_trace)) {
    show_trace = false;
  }

  if (!optional::has_value(tree_algorithm)) {
    tree_algorithm = "mst_prim";
  }

  if (!optional::has_value(seeds)) {
    seeds = std::vector<int>();
  }

  // Normalize weights if provided
  if (optional::has_value(weights) && optional::value(weights).size() > 0) {
    auto w = optional::value(weights);
    weights = w / w.sum() * w.size();
  }
}

//! @brief Validates parametric method.
inline void
FitControlsConfig::check_parametric_method() const
{
  if (optional::has_value(parametric_method)) {
    const auto& method = optional::value(parametric_method);
    if (!tools_stl::is_member(method, { "itau", "mle" })) {
      throw std::runtime_error("parametric_method should be mle or itau");
    }
  }
}

//! @brief Validates nonparametric method.
inline void
FitControlsConfig::check_nonparametric_method() const
{
  if (optional::has_value(nonparametric_method)) {
    const auto& method = optional::value(nonparametric_method);
    if (!tools_stl::is_member(method, { "constant", "linear", "quadratic" })) {
      throw std::runtime_error(
        "nonparametric_method should be constant, linear or quadratic");
    }
  }
}

//! @brief Validates nonparametric multiplier.
inline void
FitControlsConfig::check_nonparametric_mult() const
{
  if (optional::has_value(nonparametric_mult)) {
    if (optional::value(nonparametric_mult) <= 0.0) {
      throw std::runtime_error("nonparametric_mult must be positive");
    }
  }
}

//! @brief Validates nonparametric grid size.
inline void
FitControlsConfig::check_nonparametric_grid_size() const
{
  if (optional::has_value(nonparametric_grid_size)) {
    if (optional::value(nonparametric_grid_size) < 3) {
      throw std::runtime_error("nonparametric_grid_size must be at least 3");
    }
  }
}

//! @brief Validates selection criterion.
inline void
FitControlsConfig::check_selection_criterion() const
{
  if (optional::has_value(selection_criterion)) {
    const auto& criterion = optional::value(selection_criterion);
    std::vector<std::string> allowed_crits = {
      "loglik", "aic", "bic", "mbic", "mbicv"
    };
    if (!tools_stl::is_member(criterion, allowed_crits)) {
      throw std::runtime_error(
        "selection_criterion should be 'loglik', 'aic', 'bic', or 'mbic'");
    }
  }
}

//! @brief Validates psi0 parameter.
inline void
FitControlsConfig::check_psi0() const
{
  if (optional::has_value(psi0)) {
    const auto& val = optional::value(psi0);
    if ((val <= 0.0) || (val >= 1.0)) {
      throw std::runtime_error("psi0 must be in the interval (0, 1)");
    }
  }
}

//! @brief Validates tree criterion.
inline void
FitControlsConfig::check_tree_criterion() const
{
  if (optional::has_value(tree_criterion)) {
    const auto& criterion = optional::value(tree_criterion);
    if (!tools_stl::is_member(criterion,
                              { "tau", "rho", "joe", "hoeffd", "mcor" })) {
      throw std::runtime_error("tree_criterion must be one of "
                               "'tau', 'rho', 'hoeffd', 'mcor', or 'joe'");
    }
  }
}

//! @brief Validates threshold parameter.
inline void
FitControlsConfig::check_threshold() const
{
  if (optional::has_value(threshold)) {
    const auto& val = optional::value(threshold);
    if (val < 0 || val > 1) {
      throw std::runtime_error("threshold should be in [0,1]");
    }
  }
}

//! @brief Processes and validates number of threads.
inline size_t
FitControlsConfig::process_num_threads(size_t num_threads) const
{
  // zero threads means everything is done in main thread
  if (num_threads == 1)
    num_threads = 0;

  // don't use more threads than supported by the system
  size_t max_threads = std::thread::hardware_concurrency();
  num_threads = std::min(num_threads, max_threads);

  return num_threads;
}

//! @brief Conversion factory method from FitControlsBicop
//! This can only be implemented after FitControlsBicop is fully defined
// inline FitControlsConfig FitControlsConfig::from_bicop(const
// FitControlsBicop& controls);

//! @brief Conversion factory method from FitControlsVinecop
//! This can only be implemented after FitControlsVinecop is fully defined
// inline FitControlsConfig FitControlsConfig::from_vinecop(const
// FitControlsVinecop& controls);

}
