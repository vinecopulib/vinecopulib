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

//! A class for vine copulas
class Vinecop {
public:
    Vinecop(const int& d);
    
    BicopPtr get_pair_copula(int tree, int edge);
    int get_family(int tree, int edge);
    int get_rotation(int tree, int edge);
    VecXd get_parameters(int tree, int edge);
    
    MatXd construct_d_vine_matrix(const VecXd& order);
    
private:
    int d_;
    MatXd vine_matrix_;
    std::vector<BicopPtr> pair_copulas_;
};


#endif
