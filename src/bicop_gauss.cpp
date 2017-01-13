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
