// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/hana.hpp>
#include <limits>
#include <random>
#include <stdexcept>
#include <thread>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_print.hpp>

namespace defaults {

template<typename T>
struct ParameterDefault
{
  T value;
  std::function<void(T const&)> check_fn;
  std::function<std::string(const T&)> print_fn;
};

static const auto family_set = ParameterDefault<std::vector<BicopFamily>>{
  bicop_families::all /* default value */,
  [](const std::vector<BicopFamily>& v) { /* check method */ },
  [](const std::vector<BicopFamily>& v) { /* print method */ }
};

// other default parameters ...
} // end namespace defaults

#define FITCONTROLSBICOP_PARAMS                                                \
  X(family_set)                                                                \
  X(parametric_method)                                                         \
  // ....

struct FitControlsBicop
{
private:
#define X(name) decltype(defaults::name.value) name##_ = defaults::name.value;
  FITCONTROLSBICOP_PARAMS
#undef X

#define X(name)                                                                \
  const auto& name() const { return name##_; }                                 \
  FitControlsBicop& name(const decltype(defaults::name.value)& v)              \
  {                                                                            \
    if (!defaults::name.check_fn(v))                                           \
      throw std::invalid_argument("invalid value for " #name);                 \
    name##_ = v;                                                               \
    return *this;                                                              \
  }
  FITCONTROLSBICOP_PARAMS
#undef X

  std::string str() const
  {
    std::ostringstream oss;
#define X(name)                                                                \
  oss << #name << ": " << defaults::name.print_fn(name##_) << "\n";
    FITCONTROLSBICOP_PARAMS
#undef X
    return oss.str();
  }
};

namespace vinecopulib {
namespace hana = boost::hana;

// -------------- Bicop descriptors --------------
inline auto const&
bicop_fields()
{
  static const auto t = hana::make_tuple(
    // member ptr,        default-factory, checker-fn,  label, printer-tag
    hana::make_tuple(
      &FitControls::family_set,
      [] { return bicop_families::all; },
      [](auto const& fs) { FitControls::check_family_set(fs); },
      "Family set: ",
      tools_print::PrintFamilies{}),

    hana::make_tuple(
      &FitControls::parametric_method,
      [] { return std::string("mle"); },
      [](auto const& s) { FitControls::check_parametric_method(s); },
      "Parametric method: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::nonparametric_method,
      [] { return std::string("constant"); },
      [](auto const& s) { FitControls::check_nonparametric_method(s); },
      "Nonparametric method: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::nonparametric_mult,
      [] { return 1.0; },
      [](double x) { FitControls::check_nonparametric_mult(x); },
      "Nonparametric multiplier: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::nonparametric_grid_size,
      [] { return size_t{ 30 }; },
      [](size_t n) { FitControls::check_nonparametric_grid_size(n); },
      "Nonparametric grid size: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::selection_criterion,
      [] { return std::string("aic"); },
      [](auto const& s) { FitControls::check_selection_criterion(s); },
      "Selection criterion: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::weights,
      [] { return Eigen::VectorXd(); },
      [](Eigen::VectorXd const&) { /* NO_CHECK */ },
      "Weights: ",
      tools_print::PrintWeights{}),

    hana::make_tuple(
      &FitControls::preselect_families,
      [] { return true; },
      [](bool) { /* NO_CHECK */ },
      "Preselect families: ",
      tools_print::PrintYesNo{}),

    hana::make_tuple(
      &FitControls::psi0,
      [] { return 0.9; },
      [](double x) { FitControls::check_psi0(x); },
      "mBIC prior probability: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::allow_rotations,
      [] { return true; },
      [](bool) { /* NO_CHECK */ },
      "Allow rotations: ",
      tools_print::PrintYesNo{}),

    hana::make_tuple(
      &FitControls::num_threads,
      [] { return size_t{ 0 }; }, // you post-process → max(1,·)
      [](size_t) { /* NO_CHECK */ },
      "Number of threads: ",
      tools_print::PrintSkip{}));
  return t;
}

// -------------- Vine-only descriptors --------------
inline auto const&
vinecop_fields()
{
  static const auto t = hana::make_tuple(
    hana::make_tuple(
      &FitControls::trunc_lvl,
      [] { return std::numeric_limits<size_t>::max(); },
      [](size_t) { /* NO_CHECK */ },
      "Truncation level: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::tree_criterion,
      [] { return std::string("tau"); },
      [](auto const& s) { FitControls::check_tree_criterion(s); },
      "Tree criterion: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::threshold,
      [] { return 0.0; },
      [](double x) { FitControls::check_threshold(x); },
      "Threshold: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::select_trunc_lvl,
      [] { return false; },
      [](bool) { /* NO_CHECK */ },
      "Select trunc lvl: ",
      tools_print::PrintYesNo{}),

    hana::make_tuple(
      &FitControls::select_threshold,
      [] { return false; },
      [](bool) { /* NO_CHECK */ },
      "Select threshold: ",
      tools_print::PrintYesNo{}),

    hana::make_tuple(
      &FitControls::select_families,
      [] { return true; },
      [](bool) { /* NO_CHECK */ },
      "Select families: ",
      tools_print::PrintYesNo{}),

    hana::make_tuple(
      &FitControls::show_trace,
      [] { return false; },
      [](bool) { /* NO_CHECK */ },
      "Show trace: ",
      tools_print::PrintYesNo{}),

    hana::make_tuple(
      &FitControls::tree_algorithm,
      [] { return std::string("mst_prim"); },
      [](auto const& s) { FitControls::check_tree_algorithm(s); },
      "Tree algorithm: ",
      tools_print::PrintDefault{}),

    hana::make_tuple(
      &FitControls::seeds,
      [] { return std::vector<int>{}; }, // lazy generation later
      [](std::vector<int> const&) { /* NO_CHECK */ },
      "Seeds: ",
      tools_print::PrintSkip{}));
  return t;
}

//! @brief Creates default controls for bivariate copula models.
inline FitControls
FitControls::defaults_bicop()
{
  FitControls controls;
  hana::for_each(bicop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto make = hana::at_c<1>(d);
    (controls.*ptr) = make();
  });
  return controls;
}

//! @brief Creates default controls for vine copula models.
inline FitControls
FitControls::defaults_vinecop()
{
  FitControls controls = FitControls::defaults_bicop();
  hana::for_each(vinecop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto make = hana::at_c<1>(d);
    (controls.*ptr) = make();
  });
  return controls;
}

//! @brief Validates and applies defaults for unset values for Bicop
inline void
FitControls::validate_and_set_defaults_bicop()
{
  hana::for_each(bicop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto make = hana::at_c<1>(d);
    auto check = hana::at_c<2>(d);

    auto& opt = this->*ptr;
    if (!opt)
      opt = make();
    check(opt.value());
  });

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
  // bicop first
  validate_and_set_defaults_bicop();

  // vine-only overlay + checks
  hana::for_each(vinecop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto make = hana::at_c<1>(d);
    auto check = hana::at_c<2>(d);

    auto& opt = this->*ptr;
    if (!opt)
      opt = make();
    check(opt.value());
  });

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
FitControls::check_family_set(const std::vector<BicopFamily>& family_set)
{
  // Non empty set
  if (family_set.size() == 0) {
    throw std::runtime_error("family_set must not be empty");
  }
}

//! @brief Validates parametric method.
inline void
FitControls::check_parametric_method(const std::string& method)
{
  if (!tools_stl::is_member(method, { "itau", "mle" })) {
    throw std::runtime_error("parametric_method should be mle or itau");
  }
}

//! @brief Validates nonparametric method.
inline void
FitControls::check_nonparametric_method(const std::string& method)
{
  if (!tools_stl::is_member(method, { "constant", "linear", "quadratic" })) {
    throw std::runtime_error(
      "nonparametric_method should be constant, linear or quadratic");
  }
}

//! @brief Validates nonparametric multiplier.
inline void
FitControls::check_nonparametric_mult(double mult)
{
  if (mult <= 0.0) {
    throw std::runtime_error("nonparametric_mult must be positive");
  }
}

//! @brief Validates nonparametric grid size.
inline void
FitControls::check_nonparametric_grid_size(size_t grid_size)
{
  if (grid_size < 3) {
    throw std::runtime_error("nonparametric_grid_size must be at least 3");
  }
}

//! @brief Validates selection criterion.
inline void
FitControls::check_selection_criterion(const std::string& criterion)
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
FitControls::check_psi0(double psi0)
{
  if ((psi0 <= 0.0) || (psi0 >= 1.0)) {
    throw std::runtime_error("psi0 must be in the interval (0, 1)");
  }
}

//! @brief Validates tree criterion.
inline void
FitControls::check_tree_criterion(const std::string& criterion)
{
  if (!tools_stl::is_member(criterion,
                            { "tau", "rho", "joe", "hoeffd", "mcor" })) {
    throw std::runtime_error("tree_criterion must be one of "
                             "'tau', 'rho', 'hoeffd', 'mcor', or 'joe'");
  }
}

//! @brief Validates tree algorithm.
inline void
FitControls::check_tree_algorithm(const std::string& algorithm)
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
FitControls::check_threshold(double threshold)
{
  if (threshold < 0 || threshold > 1) {
    throw std::runtime_error("threshold should be in [0,1]");
  }
}

inline void
FitControls::check_weights_size(const Eigen::MatrixXd& data,
                                const Eigen::VectorXd& weights)
{
  if ((weights.size() > 0) && (weights.size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! @brief Processes and validates number of threads.
inline size_t
FitControls::process_num_threads(size_t num_threads)
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
FitControls::normalize_weights(const Eigen::VectorXd& weights)
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

  hana::for_each(bicop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto label = hana::at_c<3>(d);
    auto printer = hana::at_c<4>(d);
    tools_print::print_field(os, label, this->*ptr, printer);
  });

  tools_print::print_field(
    os, "Number of threads: ", num_threads, tools_print::PrintThreads{});

  return os.str();
}

inline std::string
FitControls::str_vinecop() const
{
  std::ostringstream os;

  hana::for_each(bicop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto label = hana::at_c<3>(d);
    auto printer = hana::at_c<4>(d);
    tools_print::print_field(os, label, this->*ptr, printer);
  });

  hana::for_each(vinecop_fields(), [&](auto const& d) {
    auto ptr = hana::at_c<0>(d);
    auto label = hana::at_c<3>(d);
    auto printer = hana::at_c<4>(d);
    tools_print::print_field(os, label, this->*ptr, printer);
  });

  tools_print::print_field(
    os, "Number of threads: ", num_threads, tools_print::PrintThreads{});

  return os.str();
}
}
