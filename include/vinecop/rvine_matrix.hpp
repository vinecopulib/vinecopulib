// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "misc/tools_eigen.hpp"

namespace vinecopulib
{
    //! @brief A class for regular vine matrices.
    //! 
    //! A regular vine (R-vine) matrix encodes the structure of a vine.
    //!     An examplary matrix is
    //! ```
    //! 1 1 1 1
    //! 2 2 2 0
    //! 3 3 0 0
    //! 4 0 0 0
    //! ```
    //! which encodes the following pair-copulas:
    //! ```
    //! | tree | edge | pair-copulas   |
    //! |------|------|----------------|
    //! | 0    | 0    | `(4, 1)`       |
    //! |      | 1    | `(3, 1)`       |
    //! |      | 2    | `(2, 1)`       |
    //! | 1    | 0    | `(4, 2; 1)`    |
    //! |      | 1    | `(3, 2; 1)`    |
    //! | 2    | 0    | `(4, 3; 2, 1)` |
    //! ```
    //! Denoting by `M[i][j]` the matrix entry in row `i` and column `j` (starting at
    //! 0), the pair-copula index for edge `e` in tree `t` of a `d` dimensional vine is
    //! `(M[d - 1 - t][e], M[t][e]; M[t - 1][e], ..., M[0][e])`. Less formally,
    //! 1. Start with the counter-diagonal element of column `e` (first conditioned 
    //!    variable).
    //! 2. Jump up to the element in row `t` (second conditioned variable).
    //! 3. Gather all entries further up in column `e` (conditioning set).
    //! 
    class RVineMatrix
    {
    public:
        RVineMatrix() {}
        RVineMatrix(const Eigen::MatrixXi& matrix);

        Eigen::MatrixXi get_matrix() const;
        Eigen::VectorXi get_order() const;
        Eigen::MatrixXi in_natural_order() const;
        Eigen::MatrixXi get_max_matrix() const;
        MatrixXb get_needed_hfunc1() const;
        MatrixXb get_needed_hfunc2() const;

        static Eigen::MatrixXi construct_d_vine_matrix(const Eigen::VectorXi& order);

    private:
        int d_;
        Eigen::MatrixXi matrix_;
    };

    int relabel_one(const int& x, const Eigen::VectorXi& old_labels, const Eigen::VectorXi& new_labels);
    Eigen::MatrixXi relabel_elements(const Eigen::MatrixXi& matrix, const Eigen::VectorXi& new_labels);
}
