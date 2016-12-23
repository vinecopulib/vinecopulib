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

#include "include/bicop_indep.hpp"

// constructor
IndepBicop::IndepBicop()
{
    family_ = 0;
    rotation_ = 0;
    parameters_ = VecXd::Zero(1);
    parameter_bounds_ = MatXd::Zero(1, 2);
}

IndepBicop::IndepBicop(const VecXd & parameters)
{
    family_ = 0;
    rotation_ = 0;
    parameters_ = parameters;
    parameter_bounds_ = MatXd::Zero(1, 2);
}

IndepBicop::IndepBicop(const VecXd& parameters, const int& rotation)
{
    family_ = 0;
    rotation_ = rotation;
    parameters_ = parameters;
    parameter_bounds_ = MatXd::Zero(1, 2);
}

// PDF
VecXd IndepBicop::pdf(const MatXd& u)
{
    return VecXd::Ones(u.rows());
}

// hfunctions: the conditioning variable is put second
VecXd IndepBicop::hfunc1(const MatXd& u)
{
    VecXd v = u.col(1);
    return v;
}

VecXd IndepBicop::hfunc2(const MatXd& u)
{
    VecXd v = u.col(0);
    return v;
}

VecXd IndepBicop::tau_to_par(const double __attribute__((unused))& tau)
{
    return VecXd::Zero(1);
}

double IndepBicop::par_to_tau(const VecXd __attribute__((unused))& parameters)
{
    return 0.0;
}

VecXd IndepBicop::hinv1(const MatXd& u)
{
    return u.col(1);
}

VecXd IndepBicop::hinv2(const MatXd& u)
{
    return u.col(0);
}
