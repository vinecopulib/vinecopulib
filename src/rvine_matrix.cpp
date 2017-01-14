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
    max_matrix_ = to_max_matrix(matrix);
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

//! Reorder R-vine matrix to natural order (d, ..., 1 on the diagonal)
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


MatXi RVineMatrix::to_max_matrix(const MatXi& matrix) 
{
    int d = matrix.rows();
    
    MatXi no_mat = to_natural_order(matrix);
    MatXi max_matrix = no_mat;
    for (int i = d - 2; i > -1; --i) {
        for (int j = 0 ; j < i + 1; ++j) {
            max_matrix(i, j) = max_matrix.block(i, j, 2, 1).maxCoeff();
        }
    }
    
    return relabel_elements(max_matrix, matrix.diagonal());
}

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
