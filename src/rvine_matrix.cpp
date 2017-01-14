/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "include/rvine_matrix.hpp"


RVineMatrix::RVineMatrix(const MatXi& matrix)
{
    d_ = matrix.rows();
    matrix_ = matrix;
    no_matrix_ = to_natural_order(matrix);
    max_matrix_ = to_max_matrix(no_matrix_);
    needed_hfunc1_ = compute_needed_hfunc1(no_matrix_);
    needed_hfunc2_ = compute_needed_hfunc2(no_matrix_);
}

MatXi RVineMatrix::get_matrix()
{
    return matrix_;
}

MatXi RVineMatrix::get_no_matrix()
{
    return no_matrix_;
}

MatXi RVineMatrix::get_max_matrix()
{
    return max_matrix_;
}

MatXb RVineMatrix::get_needed_hfunc1()
{
    return needed_hfunc1_;
}

MatXb RVineMatrix::get_needed_hfunc2()
{
    return needed_hfunc2_;
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
MatXi RVineMatrix::to_natural_order(const MatXi& matrix) 
{
    
    // create vector of new variable labels: d, ..., 1
    int d = matrix.rows();
    std::vector<int> ivec(d);    
    std::iota(std::begin(ivec), std::end(ivec), 1); // fill with 1, ..., d
    std::reverse(std::begin(ivec), std::end(ivec));
    Eigen::Map<VecXi> new_labels(&ivec[0], d);     // convert to VecXi
    
    return relabel_elements(matrix, new_labels);
}

//! Convert to maximum matrix
//! 
//! The maximum matrix is derived from an R-vine matrix by iteratively computing
//! the (elementwise) maximum of a row and the row below (starting from the 
//! bottom). It is used in estimation and evaluation algorithms to find the right
//! pseudo observations for an edge.
//! 
//! @parameter no_matrix initial R-vine matrix, assumed to be in natural order.
//! 
//! @return An Eigen::MatrixXi containing the maximum matrix
MatXi RVineMatrix::to_max_matrix(const MatXi& no_matrix) 
{
    int d = no_matrix.rows();
    MatXi max_matrix = no_matrix;
    for (int i = d - 2; i > -1; --i) {
        for (int j = 0 ; j < i + 1; ++j) {
            max_matrix(i, j) = max_matrix.block(i, j, 2, 1).maxCoeff();
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
MatXb RVineMatrix::compute_needed_hfunc1(const MatXi& no_matrix) 
{
    int d = no_matrix.rows();
    MatXb needed_hfunc1 = MatXb::Constant(d, d, false);
    needed_hfunc1.block(1, 0, d - 1, 1) = MatXb::Constant(d - 1, 1, true);
        
    MatXi max_matrix = to_max_matrix(no_matrix);
    for (int i = 1; i < d - 1; ++i) {
        int j = d - i;
        // fill row j with true below the diagonal
        needed_hfunc1.block(i, i, d - i, 1) = MatXb::Constant(d - i, 1, true);
        
        // for diagonal, check whether matrix and maximum matrix coincide
        MatXb is_mat_j = (no_matrix.block(i, 0, 1, i).array() == j);
        MatXb is_max_j = (max_matrix.block(i, 0, 1, i).array() == j);
        needed_hfunc1(i, i) = (is_mat_j.array() && is_max_j.array()).any();
    }
    
    return needed_hfunc1;
}

MatXb RVineMatrix::compute_needed_hfunc2(const MatXi& no_matrix) 
{
    int d = no_matrix.rows();
    MatXb needed_hfunc2 = MatXb::Constant(d, d, false);
    
    MatXi max_matrix = to_max_matrix(no_matrix);
    for (int i = 1; i < d - 1; ++i) {
        int j = d - i;
        // for diagonal, check whether matrix and maximum matrix coincide
        MatXb is_mat_j = (no_matrix.block(i, 0, d - i, i).array() == j);
        MatXb is_max_j = (max_matrix.block(i, 0, d - i, i).array() == j);
        MatXb is_different = (!is_mat_j.array() && is_max_j.array());
        needed_hfunc2.block(i, i, d - i, 1) = is_different.rowwise().any();
    }    
    return needed_hfunc2;
}
//! @}

// translates matrix_entry from old to new labels
int relabel_one(const int& x, const VecXi& old_labels, const VecXi& new_labels)
{
    for (unsigned int i = 0; i < old_labels.size(); ++i) {
        if (x == old_labels[i])
            return new_labels[i];
    }
    return 0;
}

// relabels all elements of the matrix (upper triangle assumed to be 0)
MatXi relabel_elements(const MatXi& matrix, const VecXi& new_labels)
{
    int d = matrix.rows();
    VecXi old_labels = matrix.diagonal();
    MatXi new_matrix = MatXi::Zero(d, d);
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < i + 1; ++j) {
            new_matrix(i, j) = relabel_one(matrix(i, j), old_labels, new_labels);
        }
    }
    
    return new_matrix;
}
