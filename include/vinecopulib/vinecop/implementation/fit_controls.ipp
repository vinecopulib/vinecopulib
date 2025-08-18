// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/random/seed_seq.hpp>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vinecopulib/misc/tools_stl.hpp>

//! @file vinecop/implementation/fit_controls.ipp
//! @brief Fit controls for Vinecop class (Implementation).

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {

//! Default constructor.
inline FitControlsVinecop::FitControlsVinecop()
  : FitControlsBicop(FitControlsConfig::vinecop_defaults())
{
}

//! @brief Construct from a FitControlsConfig.
inline FitControlsVinecop::FitControlsVinecop(FitControlsConfig cfg)
  : FitControlsBicop(std::move(cfg))
{
}

//! @brief Legacy constructor for backward compatibility.
inline FitControlsVinecop::FitControlsVinecop(
  std::vector<BicopFamily> family_set,
  std::string parametric_method,
  std::string nonparametric_method,
  double nonparametric_mult,
  size_t nonparametric_grid_size,
  size_t trunc_lvl,
  std::string tree_criterion,
  double threshold,
  std::string selection_criterion,
  const Eigen::VectorXd& weights,
  double psi0,
  bool preselect_families,
  bool select_trunc_lvl,
  bool select_threshold,
  bool select_families,
  bool show_trace,
  size_t num_threads,
  std::string tree_algorithm,
  bool allow_rotations,
  std::vector<int> seeds)
{
  auto cfg = FitControlsConfig::vinecop_defaults();
  cfg.family_set = family_set;
  cfg.parametric_method = parametric_method;
  cfg.nonparametric_method = nonparametric_method;
  cfg.nonparametric_mult = nonparametric_mult;
  cfg.nonparametric_grid_size = nonparametric_grid_size;
  cfg.selection_criterion = selection_criterion;
  cfg.weights = weights;
  cfg.psi0 = psi0;
  cfg.preselect_families = preselect_families;
  cfg.allow_rotations = allow_rotations;
  cfg.num_threads = num_threads;
  cfg.trunc_lvl = trunc_lvl;
  cfg.tree_criterion = tree_criterion;
  cfg.threshold = threshold;
  cfg.select_trunc_lvl = select_trunc_lvl;
  cfg.select_threshold = select_threshold;
  cfg.select_families = select_families;
  cfg.show_trace = show_trace;
  cfg.tree_algorithm = tree_algorithm;
  cfg.seeds = seeds;
  
  config() = cfg;
  config().validate_and_set_defaults();
}

//! @brief Legacy constructor from FitControlsBicop.
inline FitControlsVinecop::FitControlsVinecop(const FitControlsBicop& controls,
                                              size_t trunc_lvl,
                                              std::string tree_criterion,
                                              double threshold,
                                              bool select_trunc_lvl,
                                              bool select_threshold,
                                              bool select_families,
                                              bool show_trace,
                                              std::string tree_algorithm,
                                              std::vector<int> seeds)
  : FitControlsBicop(controls)
{
  auto& cfg = config();
  cfg.trunc_lvl = trunc_lvl;
  cfg.tree_criterion = tree_criterion;
  cfg.threshold = threshold;
  cfg.select_trunc_lvl = select_trunc_lvl;
  cfg.select_threshold = select_threshold;
  cfg.select_families = select_families;
  cfg.show_trace = show_trace;
  cfg.tree_algorithm = tree_algorithm;
  cfg.seeds = seeds;
  cfg.validate_and_set_defaults();
}

//! @brief Gets the truncation level (deprecated).
inline size_t
FitControlsVinecop::get_truncation_level() const
{
  return get_trunc_lvl();
}

//! @brief Gets the truncation level.
inline size_t
FitControlsVinecop::get_trunc_lvl() const
{
  return optional::value(config().trunc_lvl);
}

//! @brief Gets the tree criterion.
inline std::string
FitControlsVinecop::get_tree_criterion() const
{
  return optional::value(config().tree_criterion);
}

//! @brief Gets the threshold parameter.
inline double
FitControlsVinecop::get_threshold() const
{
  return optional::value(config().threshold);
}

//! @brief Gets whether to show a trace during fitting.
inline bool
FitControlsVinecop::get_show_trace() const
{
  return optional::value(config().show_trace);
}

//! @brief Gets whether to select the truncation level automatically (deprecated).
inline bool
FitControlsVinecop::get_select_truncation_level() const
{
  return get_select_trunc_lvl();
}

//! @brief Gets whether to select the truncation level automatically.
inline bool
FitControlsVinecop::get_select_trunc_lvl() const
{
  return optional::value(config().select_trunc_lvl);
}

//! @brief Gets whether to select the threshold automatically.
inline bool
FitControlsVinecop::get_select_threshold() const
{
  return optional::value(config().select_threshold);
}

//! @brief Gets whether to select the families automatically.
inline bool
FitControlsVinecop::get_select_families() const
{
  return optional::value(config().select_families);
}

//! @brief Gets whether sparse selection is needed.
inline bool
FitControlsVinecop::needs_sparse_select() const
{
  return (get_select_trunc_lvl() | get_select_threshold());
}

//! @brief Gets the fit controls for bivariate fitting.
inline FitControlsBicop
FitControlsVinecop::get_fit_controls_bicop() const
{
  return *static_cast<const FitControlsBicop*>(this);
}

//! @brief Gets the maximum spanning tree algorithm.
inline std::string
FitControlsVinecop::get_tree_algorithm() const
{
  return optional::value(config().tree_algorithm);
}

//! @brief Gets the random seeds for the random number generator.
inline std::vector<int>
FitControlsVinecop::get_seeds() const
{
  return optional::value(config().seeds);
}

//! @brief Gets the random number generator.
inline boost::random::mt19937
FitControlsVinecop::get_rng() const
{
  auto seeds = get_seeds();
  if (seeds.empty()) {
    rng_.seed(static_cast<boost::random::mt19937::result_type>(std::random_device{}()));
  } else {
    boost::random::seed_seq seq(seeds.begin(), seeds.end());
    rng_.seed(seq);
  }
  return rng_;
}

//! @brief Sets the truncation level (deprecated).
inline void
FitControlsVinecop::set_truncation_level(size_t trunc_lvl)
{
  set_trunc_lvl(trunc_lvl);
}

//! @brief Sets the truncation level.
inline void
FitControlsVinecop::set_trunc_lvl(size_t trunc_lvl)
{
  config().trunc_lvl = trunc_lvl;
}

//! @brief Sets the tree criterion.
inline void
FitControlsVinecop::set_tree_criterion(std::string tree_criterion)
{
  config().tree_criterion = tree_criterion;
  config().check_tree_criterion();
}

//! @brief Sets the threshold parameter.
inline void
FitControlsVinecop::set_threshold(double threshold)
{
  config().threshold = threshold;
  config().check_threshold();
}

//! @brief Sets whether to show a trace during fitting.
inline void
FitControlsVinecop::set_show_trace(bool show_trace)
{
  config().show_trace = show_trace;
}

//! @brief Sets whether to select the truncation level automatically (deprecated).
inline void
FitControlsVinecop::set_select_truncation_level(bool select_trunc_lvl)
{
  set_select_trunc_lvl(select_trunc_lvl);
}

//! @brief Sets whether to select the truncation level automatically.
inline void
FitControlsVinecop::set_select_trunc_lvl(bool select_trunc_lvl)
{
  config().select_trunc_lvl = select_trunc_lvl;
}

//! @brief Sets whether to select the threshold automatically.
inline void
FitControlsVinecop::set_select_threshold(bool select_threshold)
{
  config().select_threshold = select_threshold;
}

//! @brief Sets whether to select the families automatically.
inline void
FitControlsVinecop::set_select_families(bool select_families)
{
  config().select_families = select_families;
}

//! @brief Sets the fit controls for bivariate fitting.
inline void
FitControlsVinecop::set_fit_controls_bicop(FitControlsBicop controls)
{
  // Copy the bicop-specific settings from the provided controls
  auto& cfg = config();
  const auto& other_cfg = controls.config();
  cfg.family_set = other_cfg.family_set;
  cfg.parametric_method = other_cfg.parametric_method;
  cfg.nonparametric_method = other_cfg.nonparametric_method;
  cfg.nonparametric_mult = other_cfg.nonparametric_mult;
  cfg.nonparametric_grid_size = other_cfg.nonparametric_grid_size;
  cfg.selection_criterion = other_cfg.selection_criterion;
  cfg.weights = other_cfg.weights;
  cfg.psi0 = other_cfg.psi0;
  cfg.preselect_families = other_cfg.preselect_families;
  cfg.allow_rotations = other_cfg.allow_rotations;
  cfg.num_threads = other_cfg.num_threads;
}

//! @brief Sets the maximum spanning tree algorithm.
inline void
FitControlsVinecop::set_tree_algorithm(std::string tree_algorithm)
{
  config().tree_algorithm = tree_algorithm;
}

//! @brief Sets the random seeds for the random number generator.
inline void
FitControlsVinecop::set_seeds(std::vector<int> seeds)
{
  config().seeds = seeds;
}

//! @brief Summarizes the controls into a string (can be used for printing).
inline std::string
FitControlsVinecop::str() const
{
  std::stringstream controls_str;

  controls_str << str_internal(false);
  controls_str << "Truncation level: "
               << (get_trunc_lvl() == std::numeric_limits<size_t>::max()
                     ? "none (default)"
                     : std::to_string(get_trunc_lvl()))
               << std::endl;
  controls_str << "Tree criterion: " << get_tree_criterion() << std::endl;
  controls_str << "Threshold: " << get_threshold() << std::endl;
  controls_str << "Select truncation level: "
               << static_cast<std::string>(get_select_trunc_lvl() ? "yes"
                                                                  : "no")
               << std::endl;
  controls_str << "Select threshold: "
               << static_cast<std::string>(get_select_threshold() ? "yes"
                                                                  : "no")
               << std::endl;
  controls_str << "Select families: "
               << static_cast<std::string>(get_select_families() ? "yes" : "no")
               << std::endl;
  controls_str << "Show trace: "
               << static_cast<std::string>(get_show_trace() ? "yes" : "no")
               << std::endl;
  controls_str << "Number of threads: "
               << (get_num_threads() == 0 ? 1 : get_num_threads()) << std::endl;
  controls_str << "MST algorithm: " << get_tree_algorithm() << std::endl;
  return controls_str.str().c_str();
}

}
