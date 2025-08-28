// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <limits>
#include <string>
#include <thread>
#include <vector>
#include <vinecopulib/misc/tools_optional.hpp>

// Forward declaration - we need the full include later for implementations
namespace vinecopulib {
enum class BicopFamily;
class FitControlsBicop;
class FitControlsVinecop;
}

namespace vinecopulib {

//! @brief controls options for bivariate and vine copula models.
//!
//! @details This struct provides a flexible way to controlsure copula fitting parameters.
//! Each field is optional, and default values are applied if the field
//! is not set. Includes built-in validation when values are set.
struct FitControls
{
  //! @brief The set of copula families to consider (if empty, then all families are
  //! included).
  optional::optional<std::vector<BicopFamily>> family_set;

  //! @brief The fit method for parametric families; possible choices: `"mle"`,
  //! `"itau"`.
  optional::optional<std::string> parametric_method;

  //! @brief The fit method for the local-likelihood nonparametric family (TLLs);
  //! possible choices: `"constant"`, `"linear"`, `"quadratic"`.
  optional::optional<std::string> nonparametric_method;

  //! @brief A factor with which the smoothing parameters are multiplied.
  optional::optional<double> nonparametric_mult;

  //! @brief Grid size for nonparametric estimation.
  optional::optional<size_t> nonparametric_grid_size;

  //! @brief The selection criterion (`"loglik"`, `"aic"` or `"bic"`) for
  //! the pair copula families.
  optional::optional<std::string> selection_criterion;

  //! @brief A vector of weights for the observations.
  optional::optional<Eigen::VectorXd> weights;

  //! @brief Only for `selection_criterion = "mbic"`, prior probability of
  //! non-independence.
  optional::optional<double> psi0;

  //! @brief Whether to exclude families before fitting based on symmetry properties of
  //! the data.
  optional::optional<bool> preselect_families;

  //! @brief Whether to allow rotations for the families when doing model selection.
  optional::optional<bool> allow_rotations;

  //! @brief Number of concurrent threads to use while fitting pair copulas within a
  //! tree; never uses more than the number of concurrent threads supported by
  //! the implementation.
  optional::optional<size_t> num_threads;

  //! @brief Truncation level for truncated vines.
  optional::optional<size_t> trunc_lvl;

  //! @brief The criterion for selecting the spanning tree (`"tau"`, `"hoeffd"`,
  //! `"rho"`, and `"mcor"` implemented so far) during the tree-wise structure
  //! selection.
  optional::optional<std::string> tree_criterion;

  //! @brief For thresholded vines (0 = no threshold).
  optional::optional<double> threshold;

  //! @brief Whether the threshold parameter shall be selected automatically.
  optional::optional<bool> select_threshold;

  //! @brief Whether the truncation shall be selected automatically.
  optional::optional<bool> select_trunc_lvl;

  //! @brief Whether the families shall be selected automatically
  //!
  //! @details If true, methods will automatically select the families
  //! to be used for the pair copulas. If false, methods will simply
  //! update the parameters for the pair copulas already present in the
  //! model.
  optional::optional<bool> select_families;

  //! @brief Whether to show a trace of the building progress.
  optional::optional<bool> show_trace;

  //! @brief The algorithm for building the spanning tree (`"mst_prim"`,
  //! `"mst_kruskal"`, `"random_weighted"`, or
  //! `"random_unweighted"`) during the tree-wise structure selection.
  //! 
  //! @details `"mst_prim"` and `"mst_kruskal"` use Prim's and Kruskal's algorithms
  //! respectively to select the maximum spanning tree, maximizing
  //! the sum of the edge weights (i.e., `tree_criterion`).
  //! `"random_weighted"` and `"random_unweighted"` use Wilson's
  //! algorithm to generate a random spanning tree, either with probability
  //! proportional to the product of the edge weights (weighted) or
  //! uniformly (unweighted).
  optional::optional<std::string> tree_algorithm;

  //! @brief A vector of random seeds
  //!
  //! @details Used for random number generators
  //! in parts of the algorithms that are randomized (e.g., random tree
  //! selection).
  optional::optional<std::vector<int>> seeds;

  //! @brief Default constructor
  FitControls() = default;

  //! @brief Factory method for creating default Bicop controls.
  static FitControls defaults_bicop();

  //! @brief Factory method for creating default Vinecop controls.
  static FitControls defaults_vinecop();

  //! @brief Validates and applies defaults for unset values for Vinecop
  void validate_and_set_defaults_vinecop();

  //! @brief Validates and applies defaults for unset values for Bicop
  void validate_and_set_defaults_bicop();

  void check_family_set(const std::vector<BicopFamily>& family_set) const;
  void check_parametric_method(const std::string& method) const;
  void check_nonparametric_method(const std::string& method) const;
  void check_nonparametric_mult(double mult) const;
  void check_nonparametric_grid_size(size_t grid_size) const;
  void check_selection_criterion(const std::string& criterion) const;
  void check_psi0(double psi0) const;
  void check_tree_criterion(const std::string& criterion) const;
  void check_tree_algorithm(const std::string& algorithm) const;
  void check_threshold(double threshold) const;
  void check_weights_size(const Eigen::MatrixXd& data) const;
  //! @}

  //! @brief Processes and validates number of threads.
  size_t process_num_threads(size_t num_threads) const;

  //! @brief Normalizes weights to sum to n.
  Eigen::VectorXd normalize_weights(Eigen::VectorXd& weights) const;

  //! @brief Whether sparse selection is needed
  bool needs_sparse_select() const;

  //! @brief String representation of the controls for Bicop
  std::string str_bicop() const;

  //! @brief String representation of the controls for Vinecop
  std::string str_vinecop() const;
};

} // namespace vinecopulib

#include <vinecopulib/misc/implementation/fit_controls.ipp>