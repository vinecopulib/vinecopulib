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


#include "bicop_clayton.hpp"

// constructor
ClaytonBicop::ClaytonBicop()
{
    family_ = 3;
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Zero(1);
    parameter_bounds_ = MatXd::Zero(1, 2);
    parameter_bounds_(0, 1) = 200.0;
}

ClaytonBicop::ClaytonBicop(const VecXd& parameters)
{
    family_ = 3;
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = parameters;
    parameter_bounds_ = MatXd::Zero(1, 2);
    parameter_bounds_(0, 1) = 200.0;
}

ClaytonBicop::ClaytonBicop(const VecXd& parameters, const int& rotation)
{
    family_ = 3;
    rotation_ = rotation;
    if ((rotation == 90) | (rotation == 270))
        association_direction_ = "negative";
    else
        association_direction_ = "positive";
    parameters_ = parameters;
    parameter_bounds_ = MatXd::Zero(1, 2);
    parameter_bounds_(0, 1) = 200.0;
}

VecXd ClaytonBicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = u.array().pow(-theta);
    psi = (psi - VecXd::Ones(psi.size())) / theta;
    return psi;
}
VecXd ClaytonBicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (theta * u + VecXd::Ones(u.size()));
    psi = psi.array().pow(-1/theta);
    return psi;
}

VecXd ClaytonBicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1) * u.array().pow(-1-theta);
    return psi;
}

VecXd ClaytonBicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (1 + theta) * u.array().pow(-2-theta);
    return psi;
}

VecXd ClaytonBicop::hinv(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd hinv = u.col(0).array().pow(theta + 1.0);
    if (theta < 75) {
        hinv = u.col(1).cwiseProduct(hinv);
        hinv = hinv.array().pow(-theta/(theta + 1.0));
        VecXd x = u.col(0);
        x = x.array().pow(-theta);
        hinv = hinv - x + VecXd::Ones(x.size());
        hinv = hinv.array().pow(-1/theta);
    } else {
        hinv = hinv1_num(u);
    }
    return hinv;

}

// link between Kendall's tau and the par_bicop parameter
VecXd ClaytonBicop::tau_to_parameters(const double& tau)
{
    VecXd parameters(1);
    parameters(0) = 2 * std::fabs(tau) / (1 - std::fabs(tau));
    return parameters;
}

double ClaytonBicop::parameters_to_tau(const VecXd& parameters)
{
    double tau =  parameters(0) / (2 + std::fabs(parameters(0)));
    if ((rotation_ == 90) | (rotation_ == 270))
        tau *= -1;
    return tau;
}
