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


//! Construct a vine copula object
//! 
//! @param d tree index (starting with 1).
//! 
//! @return a d-dimensional D-vine with variable order 1, ..., d and all 
//! pair-copulas set to independence
Vinecop::Vinecop(const int& d)
{
    d_ = d;
    
    // D-vine with variable order (1, ..., d)
    VecXi order(d);
    for (int i = 0; i < d; ++i)
        order(i) = i + 1;
    RVineMatrix vine_matrix_(RVineMatrix::construct_d_vine_matrix(order));
    
    
    // all pair-copulas are independence
    pair_copulas_.reserve(d * (d - 1) / 2);
    for (unsigned int i = 0; i < pair_copulas_.size(); ++i) {
        pair_copulas_.push_back(BicopPtr(new IndepBicop));;
    }
}

//! Access to a pair copula
//! 
//! @param tree tree index (starting with 1).
//! @param edge edge index (starting with 1).
//! 
//! @return A \code std::shared_ptr (alias \code BicopPtr) to a \code Bicop 
//! object.
BicopPtr Vinecop::get_pair_copula(int tree, int edge) 
{
    if ((edge < 1) | (edge > d_ - tree)) {
        std::stringstream message;
        message << 
            "edge index out of bounds" << std::endl <<
            "allowed: 1, ..., " << d_ - tree << std::endl <<
            "actual: " << edge << std::endl << 
            "tree level: " <<  tree  << std::endl;
            throw std::runtime_error(message.str().c_str());
    }
    return pair_copulas_[(tree - 1) * d_ - tree * (tree - 1) / 2 + edge - 1];
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
