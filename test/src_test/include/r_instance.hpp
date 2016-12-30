/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecoplib.

    vinecoplib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecoplib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VINECOPLIB_R_INSTANCE_H
#define VINECOPLIB_R_INSTANCE_H

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


#endif //VINECOPLIB_R_INSTANCE_H
