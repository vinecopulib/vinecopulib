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


//! Construct a vine copula object of dimension d
//! 
//! @param d dimension of the vine copula.
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
    vine_matrix_ = RVineMatrix(RVineMatrix::construct_d_vine_matrix(order));
    
    // all pair-copulas are independence
    pair_copulas_.reserve(d * (d - 1) / 2);
    for (int i = 0; i < d * (d - 1) / 2; ++i) {
        pair_copulas_.push_back(BicopPtr(new IndepBicop));;
    }
}

//! Construct a vine copula object from a vector<BicopPtr> and structure matrix
//! 
//! @param pair_copulas pointers to Bicop objects.
//! @param matrix R-vine matrix.
//! 
//! @return a d-dimensional D-vine with variable order 1, ..., d and all 
//! pair-copulas set to independence
Vinecop::Vinecop(const std::vector<BicopPtr>& pair_copulas, const MatXi& matrix)
{
    d_ = matrix.rows();
    if ((int) pair_copulas.size() != d_ * (d_ + 1) / 2) {
        std::stringstream message;
        message << 
            "wrong size of pair_copulas; "
            "expected:" << d_ * (d_ + 1) / 2 << ", "<<
            "actual:" << pair_copulas.size() << std::endl;
            throw std::runtime_error(message.str().c_str());
    }
    // D-vine with variable order (1, ..., d)
    vine_matrix_ = RVineMatrix(matrix);
    pair_copulas_ = pair_copulas;
}

//! Access to a pair copula
//! 
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
//! 
//! @return A \code std::shared_ptr (alias \code BicopPtr) to a \code Bicop 
//! object.
BicopPtr Vinecop::get_pair_copula(int tree, int edge) 
{
    if (tree > d_ - 2) {
        std::stringstream message;
        message << 
            "tree index out of bounds" << std::endl <<
            "allowed: 0, ..., " << d_ - 2 << std::endl <<
            "actual: " << tree << std::endl;
            throw std::runtime_error(message.str().c_str());
    }
    if ((edge < 0) | (edge > d_ - tree - 2)) {
        std::stringstream message;
        message << 
            "edge index out of bounds" << std::endl <<
            "allowed: 0, ..., " << d_ - tree - 2 << std::endl <<
            "actual: " << edge << std::endl << 
            "tree level: " <<  tree  << std::endl;
            throw std::runtime_error(message.str().c_str());
    }
    return pair_copulas_[tree * d_ - tree * (tree + 1) / 2 + edge];
}

//! Get family of a pair copula
//! 
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
//! 
//! @return An \code int containing the family index.
int Vinecop::get_family(int tree, int edge) 
{
    return get_pair_copula(tree, edge)->get_family();
}

//! Get rotation of a pair copula
//! 
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
//! 
//! @return An \code int containing the rotation.
int Vinecop::get_rotation(int tree, int edge) 
{
    return get_pair_copula(tree, edge)->get_rotation();
}

//! Get parameters of a pair copula
//! 
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
//! 
//! @return An \code Eigen::VectorXd containing the parameters.
VecXd Vinecop::get_parameters(int tree, int edge) 
{
    return get_pair_copula(tree, edge)->get_parameters();
}


//! Probability density function of a vine copula
//! 
//! @param u mxd matrix of evaluation points
VecXd Vinecop::pdf(const MatXd& u)
{
    int d = u.cols();
    if (d != d_) {
        std::stringstream message;
        message << "u has wrong number of columns. " <<
                   "expected: " << d_ <<
                   ", actual: " << d << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    
    // info about the vine structure (reverse rows (!) for more natural indexing)
    VecXi order         = vine_matrix_.get_matrix().diagonal();
    MatXi no_matrix     = vine_matrix_.in_natural_order().colwise().reverse();
    MatXi max_matrix    = vine_matrix_.get_max_matrix().colwise().reverse();
    MatXb needed_hfunc1 = vine_matrix_.get_needed_hfunc1().colwise().reverse();
    MatXb needed_hfunc2 = vine_matrix_.get_needed_hfunc2().colwise().reverse();
    
    // all pair-copula densities will be multiplied at the end. initialize with
    // 1 so that unused elements have no effect
    MatXd pc_density = MatXd::Constant(d, d, 1.0);
    VecXd vine_density = VecXd::Constant(u.rows(), 1.0);  
    // temporary storage objects for h-functions (re-used for each eval point)
    MatXd direct(d, d);
    MatXd indirect(d, d);
    
    // loop through the data to keep the memory usage to a minimum; this 
    // is important for very high-dimensional settings (d >> 100)
    for (unsigned int i = 0; i < u.rows(); ++i) {
        // fill first row of direct matrix with evaluation points;
        // points have tos be reordered to correspond to natural order
        for (int j = 0; j < d; ++j)
            direct(0, j) = u(i, order(j) - 1);
    
        for (int tree = 0; tree < d - 1; ++tree) {
            for (int edge = 0; edge < d - tree - 1; ++edge) {  
                // get pair copula for this edge
                BicopPtr edge_copula = get_pair_copula(tree, edge);
                // extract evaluation point from hfunction matrices (have been
                // computed in previous tree level)
                MatXd u_e(1, 2);
                int m = max_matrix(tree, edge);
                u_e(0, 1) = direct(tree, edge);

                if (m == no_matrix(tree, edge)) {
                    u_e(0, 0) = direct(tree, d - m);
                } else {
                    u_e(0, 0) = indirect(tree, d - m);
                }

                vine_density(i) *= edge_copula->pdf(u_e)(0);
                // h-functions are only evaluated if needed in next step
                if (needed_hfunc1(tree + 1, edge))
                    direct(tree + 1, edge) = edge_copula->hfunc1(u_e)(0);
                if (needed_hfunc2(tree + 1, edge))
                    indirect(tree + 1, edge) = edge_copula->hfunc2(u_e)(0);
            }
        }
    }

    return vine_density;
}
