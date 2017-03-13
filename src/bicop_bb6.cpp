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

#include "bicop_bb6.hpp"

// constructor
Bb6Bicop::Bb6Bicop()
{
    family_ = 8;
    family_name_ = "Bb6";
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Ones(2);
    parameters_bounds_ = MatXd::Constant(2, 2, 200);
    parameters_bounds_(0, 0) = 1.0;
    parameters_bounds_(1, 0) = 1.0;
}

Bb6Bicop::Bb6Bicop(const VecXd& parameters)
{
    Bb6Bicop();
    set_parameters(parameters);
}

Bb6Bicop::Bb6Bicop(const VecXd& parameters, const int& rotation)
{
    Bb6Bicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd Bb6Bicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return std::pow((-1)*std::log(1-std::pow(1-v,theta)), delta);
    };
    return u.unaryExpr(f);
}

VecXd Bb6Bicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return 1-std::pow(1-std::exp((-1)*std::pow(v, 1/delta)), 1/theta);
    };
    return u.unaryExpr(f);
}

VecXd Bb6Bicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        double res = delta * theta *std::pow((-1)*std::log(1-std::pow(1-v,theta)),delta-1);
        return res*std::pow(1-v,theta-1)/(std::pow(1-v,theta)-1);
    };
    return u.unaryExpr(f);
}

VecXd Bb6Bicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        double tmp = std::pow(1-v,theta);
        double res = std::pow((-1)*std::log(1-tmp),delta-2)*((delta-1)*theta*tmp-(tmp+theta-1)*std::log(1-tmp));
        return res*delta*theta*std::pow(1-v, theta-2)/std::pow(tmp - 1,2);
    };
    return u.unaryExpr(f);
}

double Bb6Bicop::parameters_to_tau(const VecXd& parameters)
{
    double theta = parameters(0);
    double delta = parameters(1);
    auto f = [theta, delta](const double v) {
        return -4/(delta*theta)*std::log(1-std::pow(1-v,theta))*(1-v-std::pow(1-v,-theta)+std::pow(1-v,-theta)*v);
    };
    double tau = 1+integration_tools::integrate_zero_to_one(f);
    if ((rotation_ == 90) | (rotation_ == 270))
        tau *= -1;
    return tau;
}

