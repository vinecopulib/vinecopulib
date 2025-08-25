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

//! Configuration options for bivariate and vine copula models.
//!
//! This struct provides a flexible way to configure copula fitting parameters.
//! Each field is optional, and default values are applied if the field
//! is not set. Includes built-in validation when values are set.
struct FitControlsConfig
{
  //! The set of copula families to consider (if empty, then all families are
  //! included).
  optional::optional<std::vector<BicopFamily>> family_set;

  //! The fit method for parametric families; possible choices: `"mle"`,
  //! `"itau"`.
  optional::optional<std::string> parametric_method;

  //! The fit method for the local-likelihood nonparametric family (TLLs);
  //! possible choices: `"constant"`, `"linear"`, `"quadratic"`.
  optional::optional<std::string> nonparametric_method;

  //! A factor with which the smoothing parameters are multiplied.
  optional::optional<double> nonparametric_mult;

  //! Grid size for nonparametric estimation.
  optional::optional<size_t> nonparametric_grid_size;

  //! The selection criterion (`"loglik"`, `"aic"` or `"bic"`) for the pair
  //! copula families.
  optional::optional<std::string> selection_criterion;

  //! A vector of weights for the observations.
  optional::optional<Eigen::VectorXd> weights;

  //! Only for `selection_criterion = "mbic"`, prior probability of
  //! non-independence.
  optional::optional<double> psi0;

  //! Whether to exclude families before fitting based on symmetry properties of
  //! the data.
  optional::optional<bool> preselect_families;

  //! Whether to allow rotations for the families when doing model selection.
  optional::optional<bool> allow_rotations;

  //! Number of concurrent threads to use while fitting pair copulas within a
  //! tree; never uses more than the number of concurrent threads supported by
  //! the implementation.
  optional::optional<size_t> num_threads;

  //! Truncation level for truncated vines.
  optional::optional<size_t> trunc_lvl;

  //! The criterion for selecting the spanning tree (`"tau"`, `"hoeffd"`,
  //! `"rho"`, and `"mcor"` implemented so far) during the tree-wise structure
  //! selection.
  optional::optional<std::string> tree_criterion;

  //! For thresholded vines (0 = no threshold).
  optional::optional<double> threshold;

  //! Whether the threshold parameter shall be selected automatically.
  optional::optional<bool> select_threshold;

  //! Whether the truncation shall be selected automatically.
  optional::optional<bool> select_trunc_lvl;

  //! Whether the families shall be selected automatically, or should the method
  //! simply update the parameters for the pair copulas already present in the
  //! model.
  optional::optional<bool> select_families;

  //! Whether to show a trace of the building progress.
  optional::optional<bool> show_trace;

  //! The algorithm for building the spanning tree (`"mst_prim"`,
  //! `"mst_kruskal"`, `"random_weighted"`, or
  //! `"random_unweighted"`) during the tree-wise structure selection.
  //! `"mst_prim"` and `"mst_kruskal"` use Prim's and Kruskal's algorithms
  //! respectively to select the maximum spanning tree, maximizing
  //! the sum of the edge weights (i.e., `tree_criterion`).
  //! `"random_weighted"` and `"random_unweighted"` use Wilson's
  //! algorithm to generate a random spanning tree, either with probability
  //! proportional to the product of the edge weights (weighted) or
  //! uniformly (unweighted).
  optional::optional<std::string> tree_algorithm;

  //! A vector of random seeds for the random number generator
  //! for parts of the algorithm that are randomized (e.g., random tree
  //! selection).
  optional::optional<std::vector<int>> seeds;

  // Default constructor
  FitControlsConfig() = default;

  //! Factory methods for creating configs with appropriate defaults.
  static FitControlsConfig bicop_defaults();
  static FitControlsConfig vinecop_defaults();

  // Static factory methods from deprecated classes (defined after the full
  // types are available) static FitControlsConfig from_bicop(const
  // FitControlsBicop& controls); static FitControlsConfig from_vinecop(const
  // FitControlsVinecop& controls);

  // Validation methods - called internally when using the shim classes
  void validate_and_set_defaults();
  void check_parametric_method() const;
  void check_nonparametric_method() const;
  void check_nonparametric_mult() const;
  void check_nonparametric_grid_size() const;
  void check_selection_criterion() const;
  void check_psi0() const;
  void check_tree_criterion() const;
  void check_threshold() const;
  size_t process_num_threads(size_t num_threads) const;
};

}

#include <vinecopulib/misc/implementation/fit_controls.ipp>