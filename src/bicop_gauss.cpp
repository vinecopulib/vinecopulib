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

#include "bicop_gauss.hpp"

// constructor
GaussBicop::GaussBicop()
{
    family_ = 1;
    family_name_ = "Gaussian";
    rotation_ = 0;
    association_direction_ = "both";
    parameters_ = VecXd::Zero(1);
    parameters_bounds_ = MatXd::Ones(1, 2);
    parameters_bounds_(0, 0) = -1;
}

GaussBicop::GaussBicop(const VecXd& parameters)
{
    GaussBicop();
    set_parameters(parameters);
}

GaussBicop::GaussBicop(const VecXd& parameters, const int& rotation)
{
    GaussBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

// PDF
VecXd GaussBicop::pdf_default(const MatXd& u)
{
    // Inverse Cholesky of the correlation matrix
    double rho = double(this->parameters_(0));
    MatXd L = MatXd::Zero(2,2);
    L(0,0) = 1;
    L(1,1) = 1/sqrt(1.0-pow(rho,2.0));
    L(0,1) = -rho*L(1,1);

    // Compute copula density
    VecXd f = VecXd::Ones(u.rows());
    MatXd tmp = qnorm(u);
    f = f.cwiseQuotient(dnorm(tmp).rowwise().prod());
    tmp = tmp*L;
    f = f.cwiseProduct(dnorm(tmp).rowwise().prod());
    return f / sqrt(1.0-pow(rho,2.0));
}


// Normal h-function
VecXd GaussBicop::hfunc1_default(const MatXd& u)
{
    double rho = double(this->parameters_(0));
    VecXd h = VecXd::Zero(u.rows());
    MatXd tmp = qnorm(u);
    h = (tmp.col(1) - rho * tmp.col(0)) / sqrt(1.0 - pow(rho, 2.0));
    return pnorm(h);
}


VecXd GaussBicop::hinv1_default(const MatXd& u)
{
    double rho = double(this->parameters_(0));
    VecXd hinv = VecXd::Zero(u.rows());
    MatXd tmp = qnorm(u);
    hinv = tmp.col(1) * sqrt(1.0 - pow(rho, 2.0)) + rho * tmp.col(0);
    return pnorm(hinv);
}

VecXd GaussBicop::get_start_parameters(const double tau)
{
    return tau_to_parameters(tau);
}
