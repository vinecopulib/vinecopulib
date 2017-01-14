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

#ifndef VINECOPULIB_VINECOP_CLASS_HPP
#define VINECOPULIB_VINECOP_CLASS_HPP

#include "bicop.hpp"

typedef Eigen::MatrixXi MatXi;
typedef Eigen::VectorXi VecXi;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatXb;
// typedef Eigen::Vector<bool, Eigen::Dynamic> VecXb;

class RVineMatrix {    
public:
    RVineMatrix(const MatXi& matrix);
    MatXi get_matrix();
    MatXi get_no_matrix();
    MatXi get_max_matrix();
    MatXb get_needed_hfunc1();
    MatXb get_needed_hfunc2();

private:
    MatXi to_natural_order(const MatXi& matrix);
    MatXi to_max_matrix(const MatXi& no_matrix);
    MatXb compute_needed_hfunc1(const MatXi& no_matrix) ;
    MatXb compute_needed_hfunc2(const MatXi& no_matrix) ;

    int d_;
    MatXi matrix_;
    MatXi no_matrix_;
    MatXi max_matrix_;
    MatXb needed_hfunc1_;
    MatXb needed_hfunc2_;;
};

int relabel_one(const int& x, const VecXi& old_labels, const VecXi& new_labels);
MatXi relabel_elements(const MatXi& matrix, const VecXi& new_labels);

#endif
