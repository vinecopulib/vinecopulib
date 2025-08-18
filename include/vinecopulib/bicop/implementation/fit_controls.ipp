// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <sstream>
#include <stdexcept>
#include <thread>
#include <vinecopulib/misc/tools_stl.hpp>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {

//! @brief Default constructor.
inline FitControlsBicop::FitControlsBicop()
  : config_(FitControlsConfig::bicop_defaults())
{
  config_.validate_and_set_defaults();
}

//! @brief Construct from a FitControlsConfig.
inline FitControlsBicop::FitControlsBicop(FitControlsConfig cfg)
  : config_(std::move(cfg))
{
  config_.validate_and_set_defaults();
}

//! @brief Legacy constructor for backward compatibility.
inline FitControlsBicop::FitControlsBicop(std::vector<BicopFamily> family_set,
                                          std::string parametric_method,
                                          std::string nonparametric_method,
                                          double nonparametric_mult,
                                          size_t nonparametric_grid_size,
                                          std::string selection_criterion,
                                          const Eigen::VectorXd& weights,
                                          double psi0,
                                          bool preselect_families,
                                          bool allow_rotations,
                                          size_t num_threads)
{
  config_ = FitControlsConfig::bicop_defaults();
  config_.family_set = family_set;
  config_.parametric_method = parametric_method;
  config_.nonparametric_method = nonparametric_method;
  config_.nonparametric_mult = nonparametric_mult;
  config_.nonparametric_grid_size = nonparametric_grid_size;
  config_.selection_criterion = selection_criterion;
  config_.weights = weights;
  config_.psi0 = psi0;
  config_.preselect_families = preselect_families;
  config_.allow_rotations = allow_rotations;
  config_.num_threads = num_threads;
  config_.validate_and_set_defaults();
}

//! @brief Constructor specifying only parametric method.
inline FitControlsBicop::FitControlsBicop(std::string parametric_method)
{
  config_ = FitControlsConfig::bicop_defaults();
  config_.parametric_method = parametric_method;
  config_.validate_and_set_defaults();
}

//! @brief Constructor specifying nonparametric settings.
inline FitControlsBicop::FitControlsBicop(std::string nonparametric_method,
                                          double nonparametric_mult,
                                          size_t nonparametric_grid_size)
{
  config_ = FitControlsConfig::bicop_defaults();
  config_.nonparametric_method = nonparametric_method;
  config_.nonparametric_mult = nonparametric_mult;
  config_.nonparametric_grid_size = nonparametric_grid_size;
  config_.validate_and_set_defaults();
}

//! @brief Gets the family set.
inline std::vector<BicopFamily>
FitControlsBicop::get_family_set() const
{
  return optional::value(config_.family_set);
}

//! @brief Gets the parametric method.
inline std::string
FitControlsBicop::get_parametric_method() const
{
  return optional::value(config_.parametric_method);
}

//! @brief Gets the nonparametric method.
inline std::string
FitControlsBicop::get_nonparametric_method() const
{
  return optional::value(config_.nonparametric_method);
}

//! @brief Gets the nonparametric bandwidth multiplier.
inline double
FitControlsBicop::get_nonparametric_mult() const
{
  return optional::value(config_.nonparametric_mult);
}

//! @brief Gets the nonparametric grid size.
inline size_t
FitControlsBicop::get_nonparametric_grid_size() const
{
  return optional::value(config_.nonparametric_grid_size);
}

//! @brief Gets the selection criterion.
inline std::string
FitControlsBicop::get_selection_criterion() const
{
  return optional::value(config_.selection_criterion);
}

//! @brief Gets the observation weights.
inline Eigen::VectorXd
FitControlsBicop::get_weights() const
{
  return optional::value(config_.weights);
}

//! @brief Gets whether to preselect families.
inline bool
FitControlsBicop::get_preselect_families() const
{
  return optional::value(config_.preselect_families);
}

//! @brief Gets the baseline probability for mBIC selection.
inline double
FitControlsBicop::get_psi0() const
{
  return optional::value(config_.psi0);
}

//! @brief Gets the number of threads.
inline size_t
FitControlsBicop::get_num_threads() const
{
  return optional::value(config_.num_threads);
}

//! @brief Gets whether to allow rotations.
inline bool
FitControlsBicop::get_allow_rotations() const
{
  return optional::value(config_.allow_rotations);
}

//! @brief Sets the family set.
inline void
FitControlsBicop::set_family_set(std::vector<BicopFamily> family_set)
{
  config_.family_set = family_set;
}

//! @brief Sets the parametric method.
inline void
FitControlsBicop::set_parametric_method(std::string parametric_method)
{
  config_.parametric_method = parametric_method;
  config_.check_parametric_method();
}

//! @brief Sets the nonparametric method.
inline void
FitControlsBicop::set_nonparametric_method(std::string nonparametric_method)
{
  config_.nonparametric_method = nonparametric_method;
  config_.check_nonparametric_method();
}

//! @brief Sets the nonparametric multiplier.
inline void
FitControlsBicop::set_nonparametric_mult(double nonparametric_mult)
{
  config_.nonparametric_mult = nonparametric_mult;
  config_.check_nonparametric_mult();
}

//! @brief Sets the nonparametric grid size.
inline void
FitControlsBicop::set_nonparametric_grid_size(size_t nonparametric_grid_size)
{
  config_.nonparametric_grid_size = nonparametric_grid_size;
  config_.check_nonparametric_grid_size();
}

//! @brief Sets the selection criterion.
inline void
FitControlsBicop::set_selection_criterion(std::string selection_criterion)
{
  config_.selection_criterion = selection_criterion;
  config_.check_selection_criterion();
}

//! @brief Sets the observation weights.
inline void
FitControlsBicop::set_weights(const Eigen::VectorXd& weights)
{
  // store standardized weights (should sum up to number of observations)
  config_.weights = weights / weights.sum() * weights.size();
}

//! @brief Sets whether to preselect the families.
inline void
FitControlsBicop::set_preselect_families(bool preselect_families)
{
  config_.preselect_families = preselect_families;
}

//! @brief Sets the prior probability for mBIC.
inline void
FitControlsBicop::set_psi0(double psi0)
{
  config_.psi0 = psi0;
  config_.check_psi0();
}

//! @brief Sets the number of threads.
inline void
FitControlsBicop::set_num_threads(size_t num_threads)
{
  config_.num_threads = config_.process_num_threads(num_threads);
}

//! @brief Sets whether to allow rotations.
inline void
FitControlsBicop::set_allow_rotations(bool allow_rotations)
{
  config_.allow_rotations = allow_rotations;
}

//! @brief Summarizes the controls into a string (can be used for printing).
inline std::string
FitControlsBicop::str() const
{
  return str_internal();
}

inline std::string
FitControlsBicop::str_internal(bool print_threads) const
{
  std::stringstream controls_str;

  controls_str << "Family set: ";
  auto family_set = get_family_set();
  for (size_t j = 0; j < family_set.size(); j++) {
    if (j > 0) {
      controls_str << ", ";
    }
    controls_str << get_family_name(family_set[j]);
  }
  controls_str << std::endl;

  controls_str << "Parametric method: " << get_parametric_method() << std::endl;
  controls_str << "Nonparametric method: " << get_nonparametric_method()
               << std::endl;
  controls_str << "Nonparametric multiplier: " << get_nonparametric_mult()
               << std::endl;
  controls_str << "Nonparametric grid size: " << get_nonparametric_grid_size()
               << std::endl;
  controls_str << "Weights: "
               << static_cast<std::string>(get_weights().size() == 0 ? "no"
                                                                     : "yes")
               << std::endl;
  controls_str << "Selection criterion: " << get_selection_criterion()
               << std::endl;
  controls_str << "Preselect families: "
               << static_cast<std::string>(get_preselect_families() ? "yes"
                                                                    : "no")
               << std::endl;
  controls_str << "mBIC prior probability: " << get_psi0() << std::endl;
  if (print_threads) {
    controls_str << "Number of threads: "
                 << (get_num_threads() == 0 ? 1 : get_num_threads())
                 << std::endl;
  }
  return controls_str.str().c_str();
}

}
