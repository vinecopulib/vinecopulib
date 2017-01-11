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

#include "include/vinecop_class.hpp"



Vinecop::Vinecop(const int& d)
{
    d_ = d;
    
    // D-vine with variable order (1, ..., d)
    VecXd order(d);
    for (int i = 0; i < d; ++i)
        order(i) = i + 1;
    vine_matrix_ = construct_d_vine_matrix(order);
    
    // all pair-copulas are independence
    pair_copulas_.reserve(d * (d - 1) / 2);
    for (unsigned int i = 0; i < pair_copulas_.size(); ++i) {
        pair_copulas_.push_back(Bicop_ptr(new IndepBicop));;
    }
}

//! Access to a pair copula
//! 
//! @param tree tree index (starting with 1).
//! @param edge edge index (starting with 1).
//! 
//! @return A \code std::shared_ptr (alias \code BicopPtr) to a \code Bicop 
//! object.
Bicop_ptr Vinecop::get_pair_copula(int tree, int edge) 
{
    int pc_index = -1;
    int tree_shift = 0;
    for (int j = 1; j < d_ + 1; ++j) {
        if (tree == j) {
            pc_index = tree_shift + edge - 1;
            break;
        }
        tree_shift += d_ - j;
    }
    return pair_copulas_[pc_index];
}

int Vinecop::get_family(int tree, int edge) 
{
    return get_pair_copula(tree, edge)->get_family();
}

int Vinecop::get_rotation(int tree, int edge) 
{
    return get_pair_copula(tree, edge)->get_rotation();
}

VecXd Vinecop::get_parameters(int tree, int edge) 
{
    return get_pair_copula(tree, edge)->get_parameters();
}

//! Construct a D-vine matrix 
//! @param order order of the variables
MatXd Vinecop::construct_d_vine_matrix(const VecXd& order)
{
    int d = order.size();
    MatXd vine_matrix(d, d);
    
    for (int i = 0; i < d; ++i) {
        vine_matrix(d - i, d - i) = order(i);
    }
    
    for (int i = 0; i < d; ++i) {
        for (int j = 1; j < d - 1; ++j) {
            vine_matrix(d - i, d - j - i) = order(j);
        }
    }
    
    return vine_matrix;
}
