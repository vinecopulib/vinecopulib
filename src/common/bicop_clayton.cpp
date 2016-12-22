/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecopulib.

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


#include "include/bicop_clayton.hpp"

// constructor
ClaytonBicop::ClaytonBicop()
{
    family_ = 3;
    parameters_ = VecXd::Zero(1);
    rotation_ = 0;
}

ClaytonBicop::ClaytonBicop(double theta)
{
    family_ = 3;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = 0;
}

ClaytonBicop::ClaytonBicop(double theta, int rotation)
{
    family_ = 3;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = rotation;
}

VecXd ClaytonBicop::generator(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = u.array().pow(-theta);
    psi = (psi-VecXd::Ones(psi.size()))/theta;
    return(psi);
}
VecXd ClaytonBicop::generator_inv(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (theta*u+VecXd::Ones(u.size()));
    psi = psi.array().pow(-1/theta);
    return(psi);
}

VecXd ClaytonBicop::generator_derivative(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1)*u.array().pow(-1-theta);
    return(psi);
}

VecXd ClaytonBicop::generator_derivative2(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (1+theta)*u.array().pow(-2-theta);
    return(psi);
}

VecXd ClaytonBicop::hinv(const MatXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd hinv = u.col(0).array().pow(theta+1);
    if (theta < 75) {
        hinv = u.col(1).cwiseProduct(hinv);
        hinv = hinv.array().pow(-theta/(theta+1.0));
        VecXd x = u.col(0);
        x = x.array().pow(-theta);
        hinv = hinv-x+VecXd::Ones(x.size());
        hinv = hinv.array().pow(-1/theta);
    } else {
        hinv = hinv1_num(u);
    }
    return(hinv);

}

// link between Kendall's tau and the par_bicop parameter
double ClaytonBicop::tau_to_par(double &tau)
{
    double par = 2 * tau/(1 - std::fabs(tau));
    return(par);
}

double ClaytonBicop::par_to_tau(double &par)
{
    double tau =  par/(2 + std::fabs(par));
    return(tau);
}

MatXd ClaytonBicop::get_bounds_standard()
{
    MatXd bounds = MatXd::Zero(1,2);
    bounds(0,1) = 2e2;
    return(bounds);
}
