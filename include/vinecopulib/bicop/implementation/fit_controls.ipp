// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
#include <thread>
#include <vinecopulib/misc/tools_stl.hpp>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {

//! @brief Instantiates the controls for fitting bivariate copula models.
inline FitControlsBicop::FitControlsBicop()
: controls_(FitControls::defaults_bicop())
{
}

//! @brief Instantiates the controls for fitting bivariate copula models.
//!
//! @param family_set The set of copula families to consider (if empty, then
//!     all families are included).
//! @param parametric_method The fit method for parametric families;
//!     possible choices: `"mle"`, `"itau"`.
//! @param nonparametric_method The fit method for the local-likelihood
//!     nonparametric family (TLLs); possible choices: `"constant"`,
//!     `"linear"`, `"quadratic"`.
//! @param nonparametric_mult A factor with which the smoothing parameters
//!     are multiplied (default: 1.0).
//! @param nonparametric_grid_size The grid size for the post-estimation
//!     interpolation in nonparametric models (default: 30).
//! @param selection_criterion The selection criterion (`"loglik"`, `"aic"`
//!     or `"bic"`) for the pair copula families.
//! @param weights A vector of weights for the observations.
//! @param psi0 Only for `selection_criterion = "mbic"`, the prior probability
//!     of non-independence.
//! @param preselect_families Whether to exclude families before fitting
//!     based on symmetry properties of the data.
//! @param allow_rotations Allow rotations for the families when doing
//!     model selection (default: true).
//! @param num_threads Number of concurrent threads to use while fitting
//!     copulas for different families; never uses more than the number
//!     of concurrent threads supported by the implementation.
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
  set_family_set(family_set);
  set_parametric_method(parametric_method);
  set_nonparametric_method(nonparametric_method);
  set_nonparametric_mult(nonparametric_mult);
  set_nonparametric_grid_size(nonparametric_grid_size);
  set_selection_criterion(selection_criterion);
  set_weights(weights);
  set_preselect_families(preselect_families);
  set_allow_rotations(allow_rotations);
  set_psi0(psi0);
  set_num_threads(num_threads);
}

//! @brief Instantiates default controls except for the parameteric method.
//! @param parametric_method The fit method for parametric families;
//!     possible choices: `"mle"`, `"itau"`.
inline FitControlsBicop::FitControlsBicop(std::string parametric_method)
  : FitControlsBicop()
{
  set_parametric_method(parametric_method);
}

//! @brief Instantiates default controls except for the nonparametric method.
//! @param nonparametric_method The fit method for the local-likelihood
//!     nonparametric family (TLLs); possible choices: `"constant"`,
//!     `"linear"`, `"quadratic"`.
//! @param nonparametric_mult A factor with which the smoothing parameters
//!     are multiplied (default: 1.0).
//! @param nonparametric_grid_size The grid size for the post-estimation
//!     interpolation in nonparametric models (default: 30).
inline FitControlsBicop::FitControlsBicop(std::string nonparametric_method,
                                          double nonparametric_mult,
                                          size_t nonparametric_grid_size)
  : FitControlsBicop()
{
  set_nonparametric_method(nonparametric_method);
  set_nonparametric_mult(nonparametric_mult);
  set_nonparametric_grid_size(nonparametric_grid_size);
}

//! @brief Instantiates the controls from a uration object.
//! @param  The uration object.
inline FitControlsBicop::FitControlsBicop(const FitControls& controls)
  : controls_(controls)
{
  controls_.validate_and_set_defaults_bicop();
}

//! @name Getters and setters.
//! @{

//! @brief Gets the family set.
inline std::vector<BicopFamily>
FitControlsBicop::get_family_set() const
{
  return controls_.family_set.value();
}

//! @brief Gets the parametric method.
inline std::string
FitControlsBicop::get_parametric_method() const
{
  return controls_.parametric_method.value();
}

//! @brief Gets the nonparametric method.
inline std::string
FitControlsBicop::get_nonparametric_method() const
{
  return controls_.nonparametric_method.value();
}

//! @brief Gets the nonparametric bandwidth multiplier.
inline double
FitControlsBicop::get_nonparametric_mult() const
{
  return controls_.nonparametric_mult.value();
}

//! @brief Gets the nonparametric grid size.
inline size_t
FitControlsBicop::get_nonparametric_grid_size() const
{
  return controls_.nonparametric_grid_size.value();
}

//! @brief Gets the number of threads.
inline size_t
FitControlsBicop::get_num_threads() const
{
  return controls_.num_threads.value();
}

inline std::string
FitControlsBicop::get_selection_criterion() const
{
  return controls_.selection_criterion.value();
}

//! @brief Gets the observation weights.
inline Eigen::VectorXd
FitControlsBicop::get_weights() const
{
  return controls_.weights.value();
}

//! @brief Gets whether to preselect families.
inline bool
FitControlsBicop::get_preselect_families() const
{
  return controls_.preselect_families.value();
}

//! @brief Gets the baseline probability for mBIC selection.
inline double
FitControlsBicop::get_psi0() const
{
  return controls_.psi0.value();
}

//! @brief Gets whether to allow rotations.
inline bool
FitControlsBicop::get_allow_rotations() const
{
  return controls_.allow_rotations.value();
}

//! @brief Gets the underlying FitControls object.
inline FitControls
FitControlsBicop::get_controls() const
{
  return controls_;
}

//! @brief Sets the family set.
inline void
FitControlsBicop::set_family_set(std::vector<BicopFamily> family_set)
{
  controls_.family_set = family_set;
}

//! @brief Sets the parametric method.
inline void
FitControlsBicop::set_parametric_method(std::string parametric_method)
{
  controls_.check_parametric_method(parametric_method);
  controls_.parametric_method = parametric_method;
}

//! @brief Sets the nonparmetric method.
inline void
FitControlsBicop::set_nonparametric_method(std::string nonparametric_method)
{
  controls_.check_nonparametric_method(nonparametric_method);
  controls_.nonparametric_method = nonparametric_method;
}

//! @brief Sets the nonparametric multiplier.
inline void
FitControlsBicop::set_nonparametric_mult(double nonparametric_mult)
{
  controls_.check_nonparametric_mult(nonparametric_mult);
  controls_.nonparametric_mult = nonparametric_mult;
}

//! @brief Sets the nonparametric grid size.
inline void
FitControlsBicop::set_nonparametric_grid_size(size_t nonparametric_grid_size)
{
  controls_.check_nonparametric_grid_size(nonparametric_grid_size);
  controls_.nonparametric_grid_size = nonparametric_grid_size;
}

//! @brief Sets the selection criterion.
inline void
FitControlsBicop::set_selection_criterion(std::string selection_criterion)
{
  controls_.check_selection_criterion(selection_criterion);
  controls_.selection_criterion = selection_criterion;
}

//! @brief Sets the observation weights.
inline void
FitControlsBicop::set_weights(const Eigen::VectorXd& weights)
{
  // store standardized weights (should sum up to number of observations)
  controls_.weights = weights / weights.sum() * weights.size();
}

//! @brief Sets whether to preselect the families.
inline void
FitControlsBicop::set_preselect_families(bool preselect_families)
{
  controls_.preselect_families = preselect_families;
}

//! @brief Sets the prior probability for mBIC.
inline void
FitControlsBicop::set_psi0(double psi0)
{
  controls_.check_psi0(psi0);
  controls_.psi0 = psi0;
}

//! @brief Sets the number of threads.
inline void
FitControlsBicop::set_num_threads(size_t num_threads)
{
  controls_.num_threads = controls_.process_num_threads(num_threads);
}

//! @brief Sets whether to allow rotations.
inline void
FitControlsBicop::set_allow_rotations(bool allow_rotations)
{
  controls_.allow_rotations = allow_rotations;
}
//! @}

//! @brief Summarizes the controls into a string (can be used for printing).
inline std::string
FitControlsBicop::str() const
{
  return controls_.str_bicop();
}

}
