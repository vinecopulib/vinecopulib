// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <limits>
#include <random>
#include <stdexcept>
#include <thread>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_print.hpp>

namespace vinecopulib {

// A dummy checker that does nothing.
#define NO_CHECK(x) ((void)0)

// Bicop fields (field, default, checker, label, printer)
#define BICOP_FIELDS                                                           \
  X(family_set,                                                                \
    bicop_families::all,                                                       \
    check_family_set,                                                          \
    "Family set: ",                                                            \
    tools_print::PrintFamilies{})                                              \
  X(parametric_method,                                                         \
    "mle",                                                                     \
    check_parametric_method,                                                   \
    "Parametric method: ",                                                     \
    tools_print::PrintDefault{})                                               \
  X(nonparametric_method,                                                      \
    "constant",                                                                \
    check_nonparametric_method,                                                \
    "Nonparametric method: ",                                                  \
    tools_print::PrintDefault{})                                               \
  X(nonparametric_mult,                                                        \
    1.0,                                                                       \
    check_nonparametric_mult,                                                  \
    "Nonparametric multiplier: ",                                              \
    tools_print::PrintDefault{})                                               \
  X(nonparametric_grid_size,                                                   \
    30,                                                                        \
    check_nonparametric_grid_size,                                             \
    "Nonparametric grid size: ",                                               \
    tools_print::PrintDefault{})                                               \
  X(weights,                                                                   \
    Eigen::VectorXd(),                                                         \
    NO_CHECK,                                                                  \
    "Weights: ",                                                               \
    tools_print::PrintWeights{})                                               \
  X(selection_criterion,                                                       \
    "aic",                                                                     \
    check_selection_criterion,                                                 \
    "Selection criterion: ",                                                   \
    tools_print::PrintDefault{})                                               \
  X(preselect_families,                                                        \
    true,                                                                      \
    NO_CHECK,                                                                  \
    "Preselect families: ",                                                    \
    tools_print::PrintYesNo{})                                                 \
  X(psi0,                                                                      \
    0.9,                                                                       \
    check_psi0,                                                                \
    "mBIC prior probability: ",                                                \
    tools_print::PrintDefault{})                                               \
  X(allow_rotations,                                                           \
    true,                                                                      \
    NO_CHECK,                                                                  \
    "Allow rotations: ",                                                       \
    tools_print::PrintYesNo{})                                                 \
  X(num_threads, 0, NO_CHECK, "Number of threads: ", tools_print::PrintSkip{})

// Vinecop fields (field, default, checker, label, printer)
#define VINECOP_FIELDS                                                         \
  X(trunc_lvl,                                                                 \
    std::numeric_limits<size_t>::max(),                                        \
    NO_CHECK,                                                                  \
    "Truncation level: ",                                                      \
    tools_print::PrintDefault{})                                               \
  X(tree_criterion,                                                            \
    "tau",                                                                     \
    check_tree_criterion,                                                      \
    "Tree criterion: ",                                                        \
    tools_print::PrintDefault{})                                               \
  X(threshold,                                                                 \
    0.0,                                                                       \
    check_threshold,                                                           \
    "Threshold: ",                                                             \
    tools_print::PrintDefault{})                                               \
  X(select_trunc_lvl,                                                          \
    false,                                                                     \
    NO_CHECK,                                                                  \
    "Select trunc lvl: ",                                                      \
    tools_print::PrintYesNo{})                                                 \
  X(select_threshold,                                                          \
    false,                                                                     \
    NO_CHECK,                                                                  \
    "Select threshold: ",                                                      \
    tools_print::PrintYesNo{})                                                 \
  X(select_families,                                                           \
    true,                                                                      \
    NO_CHECK,                                                                  \
    "Select families: ",                                                       \
    tools_print::PrintYesNo{})                                                 \
  X(show_trace, false, NO_CHECK, "Show trace: ", tools_print::PrintYesNo{})    \
  X(tree_algorithm,                                                            \
    "mst_prim",                                                                \
    check_tree_algorithm,                                                      \
    "Tree algorithm: ",                                                        \
    tools_print::PrintDefault{})                                               \
  X(seeds, std::vector<int>(), NO_CHECK, "Seeds: ", tools_print::PrintSkip{})

//! @brief Creates default controls for bivariate copula models.
inline FitControls
FitControls::defaults_bicop()
{
  FitControls controls;
#define X(field, default_value, check, label, printer)                         \
  controls.field = default_value;
  BICOP_FIELDS
#undef X
  return controls;
}

//! @brief Creates default controls for vine copula models.
inline FitControls
FitControls::defaults_vinecop()
{
  FitControls controls;
#define X(field, default_value, check, label, printer)                         \
  controls.field = default_value;
  BICOP_FIELDS
  VINECOP_FIELDS
#undef X
  return controls;
}

//! @brief Validates and applies defaults for unset values for Bicop
inline void
FitControls::validate_and_set_defaults_bicop()
{
  // Start with defaults
  const auto defaults = defaults_bicop();

// Overlay
#define X(field, default_value, check, label, printer)                         \
  field = field.value_or(default_value);
  BICOP_FIELDS
#undef X

// Checks
#define X(field, default_value, check, label, printer) check(field.value());
  BICOP_FIELDS
#undef X

  // Post-processing
  num_threads = process_num_threads(num_threads.value());
  if (weights.has_value() && weights.value().size() > 0) {
    weights = normalize_weights(weights.value());
  }
}

//! @brief Validates all set values and applies defaults for unset values.
inline void
FitControls::validate_and_set_defaults_vinecop()
{
  // Overlay bicop + vine fields
#define X(field, default_value, check, label, printer)                         \
  field = field.value_or(default_value);
  BICOP_FIELDS
  VINECOP_FIELDS
#undef X

// Checks bicop + vine
#define X(field, default_value, check, label, printer) check(field.value());
  BICOP_FIELDS
  VINECOP_FIELDS
#undef X

  // Post-processing (same as bicop, plus seeds)
  num_threads = process_num_threads(num_threads.value());
  if (weights.has_value() && weights.value().size() > 0) {
    weights = normalize_weights(weights.value());
  }

  // Lazy seed generation
  if (!seeds || seeds->empty()) {
    std::random_device rd{};
    seeds.emplace(20);
    std::generate(
      seeds->begin(), seeds->end(), [&] { return static_cast<int>(rd()); });
  }
}

//! @brief Validates family set.
inline void
FitControls::check_family_set(const std::vector<BicopFamily>& family_set) const
{
  // Non empty set
  if (family_set.size() == 0) {
    throw std::runtime_error("family_set must not be empty");
  }
}

//! @brief Validates parametric method.
inline void
FitControls::check_parametric_method(const std::string& method) const
{
  if (!tools_stl::is_member(method, { "itau", "mle" })) {
    throw std::runtime_error("parametric_method should be mle or itau");
  }
}

//! @brief Validates nonparametric method.
inline void
FitControls::check_nonparametric_method(const std::string& method) const
{
  if (!tools_stl::is_member(method, { "constant", "linear", "quadratic" })) {
    throw std::runtime_error(
      "nonparametric_method should be constant, linear or quadratic");
  }
}

//! @brief Validates nonparametric multiplier.
inline void
FitControls::check_nonparametric_mult(double mult) const
{
  if (mult <= 0.0) {
    throw std::runtime_error("nonparametric_mult must be positive");
  }
}

//! @brief Validates nonparametric grid size.
inline void
FitControls::check_nonparametric_grid_size(size_t grid_size) const
{
  if (grid_size < 3) {
    throw std::runtime_error("nonparametric_grid_size must be at least 3");
  }
}

//! @brief Validates selection criterion.
inline void
FitControls::check_selection_criterion(const std::string& criterion) const
{
  std::vector<std::string> allowed_crits = {
    "loglik", "aic", "bic", "mbic", "mbicv"
  };
  if (!tools_stl::is_member(criterion, allowed_crits)) {
    throw std::runtime_error(
      "selection_criterion should be 'loglik', 'aic', 'bic', or 'mbic'");
  }
}

//! @brief Validates psi0 parameter.
inline void
FitControls::check_psi0(double psi0) const
{
  if ((psi0 <= 0.0) || (psi0 >= 1.0)) {
    throw std::runtime_error("psi0 must be in the interval (0, 1)");
  }
}

//! @brief Validates tree criterion.
inline void
FitControls::check_tree_criterion(const std::string& criterion) const
{
  if (!tools_stl::is_member(criterion,
                            { "tau", "rho", "joe", "hoeffd", "mcor" })) {
    throw std::runtime_error("tree_criterion must be one of "
                             "'tau', 'rho', 'hoeffd', 'mcor', or 'joe'");
  }
}

//! @brief Validates tree algorithm.
inline void
FitControls::check_tree_algorithm(const std::string& algorithm) const
{
  if (!tools_stl::is_member(algorithm,
                            { "mst_prim",
                              "mst_kruskal",
                              "random_weighted",
                              "random_unweighted" })) {
    throw std::runtime_error(
      "tree_algorithm must be one of 'mst_prim', 'mst_kruskal', "
      "'random_weighted', or 'random_unweighted'");
  }
}

//! @brief Validates threshold parameter.
inline void
FitControls::check_threshold(double threshold) const
{
  if (threshold < 0 || threshold > 1) {
    throw std::runtime_error("threshold should be in [0,1]");
  }
}

inline void
FitControls::check_weights_size(const Eigen::MatrixXd& data) const
{
  if ((weights.value().size() > 0) && (weights.value().size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! @brief Processes and validates number of threads.
inline size_t
FitControls::process_num_threads(size_t num_threads) const
{
  // zero threads means everything is done in main thread
  if (num_threads == 1)
    num_threads = 0;

  // don't use more threads than supported by the system
  size_t max_threads = std::thread::hardware_concurrency();
  num_threads = std::min(num_threads, max_threads);

  return num_threads;
}

//! @brief Normalizes weights to sum to n.
inline Eigen::VectorXd
FitControls::normalize_weights(Eigen::VectorXd& weights) const
{
  if (weights.size() > 0) {
    return weights / weights.sum() * weights.size();
  } else {
    throw std::runtime_error("Empty weights vector cannot be normalized.");
  }
}

inline bool
FitControls::needs_sparse_select() const
{
  return select_trunc_lvl.value_or(false) | select_threshold.value_or(false);
}

inline std::string
FitControls::str_bicop() const
{
  std::ostringstream os;

#define X(field, default_value, check, label, printer)                         \
  tools_print::print_field(os, label, field, printer);
  BICOP_FIELDS
#undef X

  tools_print::print_field(
    os, "Number of threads: ", num_threads, tools_print::PrintThreads{});

  return os.str();
}

inline std::string
FitControls::str_vinecop() const
{
  std::ostringstream os;

#define X(field, default_value, check, label, printer)                         \
  tools_print::print_field(os, label, field, printer);
  BICOP_FIELDS
  VINECOP_FIELDS
#undef X

  tools_print::print_field(
    os, "Number of threads: ", num_threads, tools_print::PrintThreads{});

  return os.str();
}

}
