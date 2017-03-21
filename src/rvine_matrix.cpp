// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "rvine_matrix.hpp"
#include "tools_stl.hpp"

namespace vinecopulib
{
    RVineMatrix::RVineMatrix(const Eigen::MatrixXi& matrix)
    {
        d_ = matrix.rows();
        // TODO: sanity checks for input matrix
        matrix_ = matrix;
    }

    Eigen::MatrixXi RVineMatrix::get_matrix() const
    {
        return matrix_;
    }

    Eigen::VectorXi RVineMatrix::get_order() const
    {
        return matrix_.colwise().reverse().diagonal().reverse();
    }

    //! Construct a D-vine matrix
    //!
    //! A D-vine is a vine where each tree is a path.
    //!
    //! @param order order of the variables
    //!
    //! @return An Eigen::MatrixXi describing the D-vine structure.
    Eigen::MatrixXi RVineMatrix::construct_d_vine_matrix(const Eigen::VectorXi& order)
    {
        int d = order.size();
        Eigen::MatrixXi vine_matrix = Eigen::MatrixXi::Constant(d, d, 0);

        for (int i = 0; i < d; ++i) {
            vine_matrix(d - 1 - i, i) = order(d - 1 - i);  // diagonal
        }

        for (int i = 1; i < d; ++i) {
            for (int j = 0; j < i; ++j) {
                vine_matrix(d - 1 - i, j) = order(i - j - 1);  // below diagonal
            }
        }

        return vine_matrix;
    }


    //! Reorder R-vine matrix to natural order
    //!
    //! Natural order means that the diagonal has entries (d, ..., 1). We convert
    //! to natural order by relabeling the variables. Most algorithms for estimation
    //! and evaluation assume that the matrix is in natural order.
    //!
    //! @parameter matrix initial R-vine matrix.
    //!
    //! @return An Eigen::MatrixXi containing the matrix in natural order.
    Eigen::MatrixXi RVineMatrix::in_natural_order() const
    {
        // create vector of new variable labels: d, ..., 1
        std::vector<int> ivec = tools_stl::seq_int(1, d_);
        tools_stl::reverse(ivec);
        Eigen::Map<Eigen::VectorXi> new_labels(&ivec[0], d_);     // convert to Eigen::VectorXi

        return relabel_elements(matrix_, new_labels);
    }

    //! Get maximum matrix
    //!
    //! The maximum matrix is derived from an R-vine matrix by iteratively computing
    //! the (elementwise) maximum of a row and the row below (starting from the
    //! bottom). It is used in estimation and evaluation algorithms to find the right
    //! pseudo observations for an edge.
    //!
    //! @parameter no_matrix initial R-vine matrix, assumed to be in natural order.
    //!
    //! @return An Eigen::MatrixXi containing the maximum matrix
    Eigen::MatrixXi RVineMatrix::get_max_matrix() const
    {
        Eigen::MatrixXi max_matrix = this->in_natural_order();
        for (int i = 0; i < d_ - 1; ++i) {
            for (int j = 0 ; j < d_ - i - 1; ++j) {
                max_matrix(i + 1, j) = max_matrix.block(i, j, 2, 1).maxCoeff();
            }
        }
        return max_matrix;
    }

    //! Obtain matrix indicating which h-functions are needed
    //!
    //! It is usually not necessary to apply both h-functions for each pair-copula.
    //!
    //! @parameter no_matrix initial R-vine matrix (assumed to be in natural order).
    //!
    //! @return An Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> indicating
    //! whether hfunc1/2 is needed for a given pair copula.
    //!
    //! @{
    MatrixXb RVineMatrix::get_needed_hfunc1() const
    {
        MatrixXb needed_hfunc1 = MatrixXb::Constant(d_, d_, false);

        Eigen::MatrixXi no_matrix = this->in_natural_order();
        Eigen::MatrixXi max_matrix = this->get_max_matrix();
        for (int i = 1; i < d_ - 1; ++i) {
            int j = d_ - i;
            MatrixXb isnt_mat_j = (no_matrix.block(0, 0, j, i).array() != j);
            MatrixXb is_max_j = (max_matrix.block(0, 0, j, i).array() == j);
            MatrixXb is_different = (isnt_mat_j.array() && is_max_j.array());
            needed_hfunc1.block(0, i, j, 1) = is_different.rowwise().any();
        }
        return needed_hfunc1;
    }

    MatrixXb RVineMatrix::get_needed_hfunc2() const
    {
        MatrixXb needed_hfunc2 = MatrixXb::Constant(d_, d_, false);
        needed_hfunc2.block(0, 0, d_ - 1, 1) = MatrixXb::Constant(d_ - 1, 1, true);
        Eigen::MatrixXi no_matrix = this->in_natural_order();
        Eigen::MatrixXi max_matrix = this->get_max_matrix();
        for (int i = 1; i < d_ - 1; ++i) {
            int j = d_ - i;
            // fill column i with true above the diagonal
            needed_hfunc2.block(0, i, d_ - i, 1) = MatrixXb::Constant(d_ - i, 1, true);
            // for diagonal, check whether matrix and maximum matrix coincide
            MatrixXb is_mat_j = (no_matrix.block(j - 1, 0, 1, i).array() == j);
            MatrixXb is_max_j = (max_matrix.block(j - 1, 0, 1, i).array() == j);
            needed_hfunc2(j - 1, i) = (is_mat_j.array() && is_max_j.array()).any();
        }

        return needed_hfunc2;
    }
    //! @}

// translates matrix_entry from old to new labels
    int relabel_one(const int& x, const Eigen::VectorXi& old_labels, const Eigen::VectorXi& new_labels)
    {
        for (unsigned int i = 0; i < old_labels.size(); ++i) {
            if (x == old_labels[i]) {
                return new_labels[i];
            }
        }
        return 0;
    }

// relabels all elements of the matrix (upper triangle assumed to be 0)
    Eigen::MatrixXi relabel_elements(const Eigen::MatrixXi& matrix, const Eigen::VectorXi& new_labels)
    {
        int d = matrix.rows();
        Eigen::VectorXi old_labels = matrix.colwise().reverse().diagonal();
        Eigen::MatrixXi new_matrix = Eigen::MatrixXi::Zero(d, d);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d - i; ++j) {
                new_matrix(i, j) = relabel_one(matrix(i, j), old_labels, new_labels);
            }
        }

        return new_matrix;
    }
}
