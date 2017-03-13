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

#include "include/bicop_parametric.hpp"
#include <RInside.h>
#include <RcppEigen.h>
#include <string>
#include <algorithm>

// A singleton class to hold the R API
class RInstance
{
    RInside R;

public:
    // family-related getters and setters
    void set_family(int family);
    void set_parameters(VecXd parameters);

    void set_rotation(int rotation);
    VecXd get_parameters();
    int get_rotation();
    int get_family();
    int get_n();

    // other methods
    RInstance();
    RInside get_R();
    VecXd eval_in_R(std::string eval_fct, int start); // evaluate a function in R
    void change_n(int n);
    void set_tau(double tau);
    double get_tau();
    MatXd get_U();

protected:
    MatXd U_;
    int n_;
    int family_;
    int rotation_;
    double tau_;
    VecXd parameters_;
};


