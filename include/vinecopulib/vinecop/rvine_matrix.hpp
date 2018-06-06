// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_stl.hpp>
#include <iostream>

namespace vinecopulib {

template<class T>
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


template<class T>
RVineMatrix<T>::RVineMatrix(size_t d) : RVineMatrix(d, d - 1) {}

template<class T>
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

template<class T>
T& RVineMatrix<T>::operator()(size_t tree, size_t edge)
{
    assert(tree < trunc_lvl_);
    assert(edge < d_ - 1 - tree);
    return mat_[edge][tree];
}

template<class T>
T RVineMatrix<T>::operator()(size_t tree, size_t edge) const 
{
    assert(tree < trunc_lvl_);
    assert(edge < d_ - 1 - tree);
    return mat_[edge][tree];
}

template<class T>
std::vector<T>& RVineMatrix<T>::operator[](size_t column) 
{
    assert(column < d_ - 1);
    return mat_[column];
}

template<class T>
std::vector<T> RVineMatrix<T>::operator[](size_t column) const 
{
    assert(column < d_ - 1);
    return mat_[column];
}

template<class T>
size_t RVineMatrix<T>::get_trunc_lvl() const 
{
    return trunc_lvl_;
}

template<class T>
size_t RVineMatrix<T>::get_dim() const
{
    return d_;
}

template<class T>
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

}

template<class T>
std::ostream& operator<<(std::ostream& os, const vinecopulib::RVineMatrix<T>& rvm) 
{  
    os << rvm.str();
    return os;  
}  
