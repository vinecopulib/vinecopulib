// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

template<class T>
class RVineMatrix {
public:
    RVineMatrix() {}
    RVineMatrix(size_t d) : d_(d)
    {
        mat_ = std::vector<std::vector<T>>(d - 1);
        for (size_t i = 0; i < mat_.size(); i++)
            mat_[i] = std::vector<T>(mat_.size() - i, 0);
    }
    
    T& operator()(size_t tree, size_t edge) {return mat_[edge][tree];}
    T operator()(size_t tree, size_t edge) const {return mat_[edge][tree];}
    //T& operator[](size_t column) {return mat_[column];}
    //T operator[](size_t column) const {return mat_[column];}

    std::string str()
    {
        std::stringstream str;
        for (size_t i = 0; i < mat_.size(); i++) {
            for (size_t j = 0; j < mat_.size() - i; j++) {
                str << (*this)(i, j) << " ";
            }
            str << std::endl;
        }
        return str.str();
    }
    
    size_t dim() const {return d_;}

private:
    size_t d_;
    std::vector< std::vector<T> > mat_;
};

}
