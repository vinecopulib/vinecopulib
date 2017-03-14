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

#include "bicop_bb8.hpp"

// constructor
Bb8Bicop::Bb8Bicop()
{
    family_ = 10;
    family_name_ = "Bb8";
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Ones(2);
    parameters_bounds_ = MatXd::Constant(2, 2, 1);
    parameters_bounds_(0, 1) = 200.0;
    parameters_bounds_(1, 0) = 0.0;
}

Bb8Bicop::Bb8Bicop(const VecXd& parameters)
{
    Bb8Bicop();
    set_parameters(parameters);
}

Bb8Bicop::Bb8Bicop(const VecXd& parameters, const int& rotation)
{
    Bb8Bicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd Bb8Bicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return -std::log((1-std::pow(1-delta*v,theta))/(1-std::pow(1-delta,theta)));
    };
    return u.unaryExpr(f);
}

VecXd Bb8Bicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return (1-std::pow(1+std::exp(-v)*(std::pow(1-delta,theta)-1),1/theta))/delta;
    };
    return u.unaryExpr(f);
}

VecXd Bb8Bicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return -delta*theta*std::pow(1-delta*v,theta-1)/(1-std::pow(1-delta*v,theta));
    };
    return u.unaryExpr(f);
}

VecXd Bb8Bicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        double tmp = std::pow(1-delta*v,theta);
        return std::pow(delta,2)*theta*std::pow(1-delta*v,theta-2)*(theta-1+tmp)/std::pow(tmp-1,2);
    };
    return u.unaryExpr(f);
}

double Bb8Bicop::parameters_to_tau(const VecXd& parameters)
{
    double theta = parameters(0);
    double delta = parameters(1);
    auto f = [theta, delta](const double t) {
        double tmp = std::pow(1-t*delta,theta);
        return std::log((tmp-1)/(std::pow(1-delta,theta)-1))*(1-t*delta-std::pow(1-t*delta,1-theta));
    };
    double tau = 1-4/(delta*theta)*integration_tools::integrate_zero_to_one(f);
    if ((rotation_ == 90) | (rotation_ == 270))
        tau *= -1;
    return tau;
}