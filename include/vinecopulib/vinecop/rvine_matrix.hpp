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
    RVineMatrix() {}
    RVineMatrix(size_t d, size_t trunc_lvl) : d_(d), trunc_lvl_(trunc_lvl)
    {
        if (d < 2)
            throw std::runtime_error("d should be greater than 1");
        //if (trunc_lvl < 1)
        //    throw std::runtime_error("trunc_lvl should be greater than 0.");
        if (trunc_lvl > d - 1)
            trunc_lvl_ = d - 1;

        mat_ = std::vector<std::vector<T>>(d - 1);
        for (size_t i = 0; i < d - 1; i++)
            mat_[i] = std::vector<T>(std::min(d - i - 1, trunc_lvl));
    }

    RVineMatrix(size_t d) : RVineMatrix(d, d - 1) {}

    T& operator()(size_t tree, size_t edge) {return mat_[edge][tree];}
    T operator()(size_t tree, size_t edge) const {return mat_[edge][tree];}
    //T& operator[](size_t column) {return mat_[column];}
    //T operator[](size_t column) const {return mat_[column];}

    std::string str() const
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

    size_t get_trunc_lvl() const {return trunc_lvl_;}
    size_t get_dim() const {return d_;}

private:
    size_t d_;
    size_t trunc_lvl_;
    std::vector<std::vector<T>> mat_;
};

}
