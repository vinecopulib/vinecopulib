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
