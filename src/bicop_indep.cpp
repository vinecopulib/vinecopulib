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

#include "bicop_indep.hpp"

// constructor
IndepBicop::IndepBicop()
{
    family_ = 0;
    family_name_ = "Independence";
    rotation_ = 0;
    association_direction_ = "none";
}

IndepBicop::IndepBicop(const VecXd& parameters)
{
    IndepBicop();
    set_parameters(parameters);
}

IndepBicop::IndepBicop(const VecXd& parameters, const int& rotation)
{
    IndepBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

// PDF
VecXd IndepBicop::pdf_default(const MatXd& u)
{
    return VecXd::Ones(u.rows());
}

// hfunctions: the conditioning variable is put second
VecXd IndepBicop::hfunc1_default(const MatXd& u)
{
    return u.col(1);
}

VecXd IndepBicop::hfunc2_default(const MatXd& u)
{
    return u.col(0);
}

VecXd IndepBicop::hinv1_default(const MatXd& u)
{
    return u.col(1);
}

VecXd IndepBicop::hinv2_default(const MatXd& u)
{
    return u.col(0);
}

VecXd IndepBicop::tau_to_parameters(const double &)
{
    VecXd pars;
    return pars;
}

double IndepBicop::parameters_to_tau(const VecXd &)
{
    return 0.0;
}

VecXd IndepBicop::get_start_parameters(const double tau)
{
    return tau_to_parameters(tau);
}

void IndepBicop::flip()
{
    // nothing to do because independence copula is radially syemmetric
}
