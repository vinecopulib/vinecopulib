// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_bb1.hpp"


// constructor
Bb1Bicop::Bb1Bicop()
{
    family_ = 7;
    family_name_ = "Bb1";
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Zero(2);
    parameters_(1) = 1;
    parameters_bounds_ = MatXd::Constant(2, 2, 200);
    parameters_bounds_(0, 0) = 0.0;
    parameters_bounds_(1, 0) = 1.0;
}

Bb1Bicop::Bb1Bicop(const VecXd& parameters)
{
    Bb1Bicop();
    set_parameters(parameters);
}

Bb1Bicop::Bb1Bicop(const VecXd& parameters, const int& rotation)
{
    Bb1Bicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd Bb1Bicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return std::pow(std::pow(v, -theta) - 1, delta);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return std::pow(std::pow(v, 1/delta) + 1, -1/theta);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return -delta * theta * std::pow(v, -(1 + theta))*std::pow(std::pow(v, -theta) - 1, delta - 1);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        double res = delta * theta * std::pow(std::pow(v, -theta) - 1, delta) / std::pow(std::pow(v, theta) - 1, 2);
        return res * (1 + delta * theta - (1 + theta) * std::pow(v, theta)) / std::pow(v, 2);
    };
    return u.unaryExpr(f);
}

double Bb1Bicop::parameters_to_tau(const VecXd& parameters)
{
    double tau = 1-2/(parameters(1) * (parameters(0) + 2));
    if ((rotation_ == 90) | (rotation_ == 270))
        tau *= -1;
    return tau;
}


/*// PDF
VecXd Bb1Bicop::pdf_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));

    VecXd f = VecXd::Ones(u.rows());
    if (theta > 1e-6 || delta > 1+1e-6)
    {
        VecXd t1 = generator(u.col(0));
        VecXd t2 = generator(u.col(1));
        VecXd ones = f;
        f = t1.cwiseProduct(t2);
        t1 = t1 + t2;
        t2 = t1.array().pow(1/delta);
        f = f.cwiseProduct((delta - 1)*theta*ones+(1+delta*theta)*t2);

        t2 = t2 + ones;
        t1 = t1.array().pow(-2+1/delta);
        t2 = t2.array().pow(-2-1/theta);
        f = f.cwiseProduct(t1).cwiseProduct(t2);

        t1 = u.col(0);
        t2 = u.col(1);
        t1 = t1.array().pow(theta);
        t2 = t2.array().pow(theta);
        t1 = t1 - ones;
        t2 = t2 - ones;
        f = f.cwiseQuotient(u.rowwise().prod()).cwiseQuotient(t1).cwiseQuotient(t2);
    }

    return f;
}

// hfunction
VecXd Bb1Bicop::hfunc1_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    VecXd t1 = generator(u.col(0));
    VecXd t2 = generator(u.col(1));
    VecXd ones = VecXd::Ones(u.rows());
    VecXd f = ones;
    t1 = t1 + t2;
    t2 = t1.array().pow(1/delta);

    t2 = t2 + ones;
    t1 = t1.array().pow(-1+1/delta);
    t2 = t2.array().pow(-1-1/theta);
    f = f.cwiseProduct(t1).cwiseProduct(t2);


    t1 = u.col(0);
    t2 = t1.unaryExpr([theta, delta](const double v){ return std::pow(std::pow(v,-theta)-1,delta-1);});
    t1 = t1.array().pow(-1-theta);
    f = f.cwiseProduct(t1).cwiseProduct(t2);
    return f;
}*/
