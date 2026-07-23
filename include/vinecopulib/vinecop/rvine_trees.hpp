// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <functional>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

//! @brief A list-of-trees decomposition of an R-vine structure.
//!
//! @details The vine is stored as one edge list per tree. Each edge holds its
//! conditioned pair `(a, b)`, its conditioning set `C`, and (optionally) the
//! fitted pair-copula. The class converts to and from the `(order,
//! struct_array)` representation used by `RVineStructure`, validating the
//! proximity condition on the way back. It is the shared primitive behind
//! `Vinecop::reorient()` and the structure finalisation during selection.
//!
//! The orientation convention is that an edge's pair-copula has its first
//! argument aligned with the conditioned variable `a`; a copula is therefore
//! flipped iff the variable placed on the matrix diagonal differs from `a`.
class RVineTrees
{
public:
  //! @brief A single edge: conditioned pair `(a, b)`, conditioning set `C`,
  //! and the associated pair-copula (independence by default).
  struct Edge
  {
    size_t a{ 0 };        //!< first conditioned variable (pair-copula arg1)
    size_t b{ 0 };        //!< second conditioned variable
    std::set<size_t> C{}; //!< conditioning set
    Bicop pair_copula{};  //!< pair-copula (independence by default)

    Edge() = default;
    Edge(size_t a_arg, size_t b_arg, std::set<size_t> c_arg)
      : a(a_arg)
      , b(b_arg)
      , C(std::move(c_arg))
    {
    }
    Edge(size_t a_arg, size_t b_arg, std::set<size_t> c_arg, Bicop pc)
      : a(a_arg)
      , b(b_arg)
      , C(std::move(c_arg))
      , pair_copula(std::move(pc))
    {
    }
  };
  using Tree = std::vector<Edge>;

  //! @brief Chooses which variable sits on the diagonal of a given column.
  //!
  //! @details Called with the column index and, for each leaf edge of the
  //! current top tree (in iteration order), the variables that may go on the
  //! diagonal (its leaf endpoints, 1 or 2 of them). Must return one of them.
  using DiagonalPolicy =
    std::function<size_t(size_t col,
                         const std::vector<std::vector<size_t>>& leaf_edges)>;

  //! @brief The result of a copula-aware `to_struct_array()`.
  struct Decomposition
  {
    std::vector<size_t> order;
    TriangularArray<size_t> struct_array;
    std::vector<std::vector<Bicop>> pair_copulas; //!< indexed `[tree][edge]`
  };

  //! @brief The augmented (line-graph) view of a tree used during conversion.
  struct AugmentedTree
  {
    std::vector<std::tuple<size_t, size_t, Edge>> edges;
    std::map<size_t, size_t> degrees;
  };

  RVineTrees()
    : d_(0)
    , trunc_lvl_(0)
    , trees_()
  {
  }

  explicit RVineTrees(const std::vector<size_t>& order,
                      const TriangularArray<size_t>& struct_array);
  RVineTrees(const std::vector<size_t>& order,
             const TriangularArray<size_t>& struct_array,
             const std::vector<std::vector<Bicop>>& pair_copulas);
  RVineTrees(size_t d, std::vector<Tree> trees);

  const std::vector<Tree>& get_trees() const { return trees_; }

  std::tuple<std::vector<size_t>, TriangularArray<size_t>> to_struct_array()
    const;
  Decomposition to_struct_array(const DiagonalPolicy& diagonal_policy) const;

  size_t get_dim() const { return d_; }
  size_t get_trunc_lvl() const { return trunc_lvl_; }

private:
  size_t d_;
  size_t trunc_lvl_;
  std::vector<Tree> trees_;

  std::map<std::pair<size_t, std::set<size_t>>, size_t> build_lookup(
    const std::vector<std::tuple<size_t, size_t, Edge>>& edges) const;
  void check_missing_vars(const std::vector<Edge>& edges, size_t d) const;
  std::vector<AugmentedTree> trees_to_augmented() const;
  Decomposition peel(const DiagonalPolicy& diagonal_policy,
                     bool carry_copulas) const;
};

}

#include <vinecopulib/vinecop/implementation/rvine_trees.ipp>
