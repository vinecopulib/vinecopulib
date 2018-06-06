// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_stl.hpp>
#include <iostream>

namespace vinecopulib {

//! @brief Right trapezoidal matrices.
//!
//! A right trapezoid is a trapezoid having two right angles. The data structure
//! `RightTrapezoid` behaves like a matrix with the structure 
//! ```
//! x x x x x
//! x x x x
//! x x x
//! x x
//! x
//! ```
//! and all other elements omitted. Such structures appear naturally in the
//! representation of a vine copula model and related algorithms. Each row 
//! corresponds to one tree in the vine, starting from the top. In each row 
//! (= tree), each column represents an edge in this tree.
//! 
//! In truncated vine models the last few rows are omitted. For example, a 
//! 3-truncated vine requires only the elements
//! ```
//! x x x x x
//! x x x x
//! x x x
//! ```
//! 
//! Only the elements indicated by `x`s are stored and can be accessed.
//! 
//! The data structure is templated and any type or class can be used to fill
//! the entries (`x`s) of the right trapezoid.
template<typename T>
class RVineMatrix {
public:
    RVineMatrix() = default;
    RVineMatrix(size_t d);
    RVineMatrix(size_t d, size_t trunc_lvl);

    T& operator()(size_t tree, size_t edge);
    T operator()(size_t tree, size_t edge) const;
    std::vector<T>& operator[](size_t column);
    std::vector<T> operator[](size_t column) const;

    size_t get_trunc_lvl() const;
    size_t get_dim() const;
    
    std::string str() const;

private:
    size_t d_;
    size_t trunc_lvl_;
    std::vector<std::vector<T>> mat_;
};


//! construct an right trapezoid of dimension `d` (the matrix has `d-1` columns 
//! and `d-1` rows).
//! @param d the dimension of the underlying vine.
template<typename T>
RVineMatrix<T>::RVineMatrix(size_t d) : RVineMatrix(d, d - 1) {}

//! construct a truncated right trapezoid (the matrix has `d-1` columns and
//! `min(trunv_lvl, d-1)` rows).
//! @param d the dimension of the vine.
//! @param trunc_lvl the truncation level.
template<typename T>
RVineMatrix<T>::RVineMatrix(size_t d, size_t trunc_lvl) : 
    d_(d), 
    trunc_lvl_(std::min(d - 1, trunc_lvl))
{
    if (d < 2)
        throw std::runtime_error("d should be greater than 1");

    mat_ = std::vector<std::vector<T>>(d - 1);
    for (size_t i = 0; i < d - 1; i++)
        mat_[i] = std::vector<T>(std::min(d - i - 1, trunc_lvl));
}

//! access one element of the trapezoid (writable).
//! @param tree the tree level.
//! @param edge the edge in this tree.
template<typename T>
T& RVineMatrix<T>::operator()(size_t tree, size_t edge)
{
    assert(tree < trunc_lvl_);
    assert(edge < d_ - 1 - tree);
    return mat_[edge][tree];
}

//! access one element of the trapezoid (non-writable).
//! @param tree the tree level.
//! @param edge the edge in this tree.
template<typename T>
T RVineMatrix<T>::operator()(size_t tree, size_t edge) const 
{
    assert(tree < trunc_lvl_);
    assert(edge < d_ - 1 - tree);
    return mat_[edge][tree];
}

//! access one column of the trapezoid (writable).
//! @param column which column to extract.
template<typename T>
std::vector<T>& RVineMatrix<T>::operator[](size_t column) 
{
    assert(column < d_ - 1);
    return mat_[column];
}

//! access one column of the trapezoid (non-writable).
//! @param column which column to extract.
template<typename T>
std::vector<T> RVineMatrix<T>::operator[](size_t column) const 
{
    assert(column < d_ - 1);
    return mat_[column];
}

//! get the truncation level of the underlying vine.
template<typename T>
size_t RVineMatrix<T>::get_trunc_lvl() const 
{
    return trunc_lvl_;
}

//! get the dimension of the underlying vine (the matrix has `d-1` columns and
//! `min(trunv_lvl, d-1)` rows).
template<typename T>
size_t RVineMatrix<T>::get_dim() const
{
    return d_;
}

//! represent RightTrapezoid as a string.
template<typename T>
std::string RVineMatrix<T>::str() const
{
    std::stringstream str;
    for (size_t i = 0; i < std::min(d_ - 1, trunc_lvl_); i++) {
        for (size_t j = 0; j < d_ - i - 1; j++) {
            str << (*this)(i, j) << " ";
        }
        str << std::endl;
    }
    return str.str();
}

} // end of namespace vinecopulib!

//! ostream method for RightTrapezoid, to be used with `std::cout`
//! @param os an output stream.
//! @param rvm an right trapezoid.
template<typename T>
std::ostream& operator<<(std::ostream& os, const vinecopulib::RVineMatrix<T>& rvm) 
{  
    os << rvm.str();
    return os;  
}  
