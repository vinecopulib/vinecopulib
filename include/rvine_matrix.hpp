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


class RVineMatrix {    
public:
    RVineMatrix(const MatXi& matrix);
    MatXi get_matrix();
    MatXi get_no_matrix();

    static MatXi to_natural_order(const MatXi& matrix);
     
private:
    int d_;
    MatXi matrix_;
    MatXi no_matrix_;
};

int relabel(const int& matrix_entry,
    const VecXi& old_labels, 
    const VecXi& new_labels);

#endif
