/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include "bicop.hpp"

typedef Eigen::MatrixXi MatXi;
typedef Eigen::VectorXi VecXi;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatXb;

class RVineMatrix {    
public:
    //! \devgroup constructors Constructors
    //! @{
    RVineMatrix() {}
    RVineMatrix(const MatXi& matrix);
    //! @}
    
    //! Getters
    //! @{
    MatXi get_matrix();
    //! @}
    
    VecXi get_order();
    MatXi in_natural_order();
    MatXi get_max_matrix();
    MatXb get_needed_hfunc1();
    MatXb get_needed_hfunc2();

    static MatXi construct_d_vine_matrix(const VecXi& order);

private:
    int d_;
    MatXi matrix_;
};

int relabel_one(const int& x, const VecXi& old_labels, const VecXi& new_labels);
MatXi relabel_elements(const MatXi& matrix, const VecXi& new_labels);

