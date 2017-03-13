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


#include "bicop_clayton.hpp"

// constructor
ClaytonBicop::ClaytonBicop()
{
    family_ = 3;
    family_name_ = "Clayton";
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Zero(1);
    parameters_bounds_ = MatXd::Zero(1, 2);
    parameters_bounds_(0, 1) = 200.0;
}

ClaytonBicop::ClaytonBicop(const VecXd& parameters)
{
    ClaytonBicop();
    set_parameters(parameters);
}

ClaytonBicop::ClaytonBicop(const VecXd& parameters, const int& rotation)
{
    ClaytonBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd ClaytonBicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (std::pow(v, -theta)-1)/theta;
    };
    return u.unaryExpr(f);
}
VecXd ClaytonBicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return std::pow(1+theta*v, -1/theta);
    };
    return u.unaryExpr(f);
}

VecXd ClaytonBicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (-1)*std::pow(v, -1-theta);
    };
    return u.unaryExpr(f);
}

VecXd ClaytonBicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (1+theta)*std::pow(v, -2-theta);
    };
    return u.unaryExpr(f);
}

// inverse h-function
VecXd ClaytonBicop::hinv1_default(const MatXd& u)
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

VecXd ClaytonBicop::get_start_parameters(const double tau)
{
    return tau_to_parameters(tau);
}

/*// PDF
VecXd ClaytonBicop::pdf_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd t1 = generator(u.col(0));
    VecXd t2 = generator(u.col(1));
    VecXd t = t1+t2;
    VecXd f = generator_inv(t);

    t1 = u.col(0).array().pow(theta);
    t2 = u.col(1).array().pow(theta);
    t = t1 + t2 - t1.cwiseProduct(t2);
    t = t.array().square();

    t1 = t1.array().pow((theta-1)/theta);
    t2 = t2.array().pow((theta-1)/theta);

    f = (1+theta)*f.cwiseQuotient(t).cwiseProduct(t1).cwiseProduct(t2);
    return f;
}

// hfunction
VecXd ClaytonBicop::hfunc1_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd t1 = generator(u.col(0));
    VecXd t2 = generator(u.col(1));
    VecXd t = t1+t2;
    VecXd f = generator_inv(t);
    f = f.array().pow(1+theta);

    t2 = u.col(0).array().pow(-1-theta);
    f = f.cwiseProduct(t2);
    return f;
}*/
