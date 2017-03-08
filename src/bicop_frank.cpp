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

#include "bicop_frank.hpp"


// constructor
FrankBicop::FrankBicop()
{
    family_ = 5;
    family_name_ = "Frank";
    rotation_ = 0;
    association_direction_ = "both";
    parameters_ = VecXd::Zero(1);
    parameters_bounds_ = MatXd::Zero(1, 2);
    parameters_bounds_(0, 0) = -200.0;
    parameters_bounds_(0, 1) = 200.0;
}

FrankBicop::FrankBicop(const VecXd& parameters)
{
    FrankBicop();
    set_parameters(parameters);
}

FrankBicop::FrankBicop(const VecXd& parameters, const int& rotation)
{
    FrankBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd FrankBicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (-1)*std::log((std::exp(-theta*v)-1)/(std::exp(-theta)-1));
    };
    return u.unaryExpr(f);
}
VecXd FrankBicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (-1/theta)*std::log(1+std::exp(-theta*v)-std::exp(-v));
    };
    return u.unaryExpr(f);
}

VecXd FrankBicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return theta*1/(1-std::exp(theta*v));
    };
    return u.unaryExpr(f);
}

VecXd FrankBicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return std::pow(theta,2)*std::exp(theta*v)/std::pow(std::exp(theta*v) - 1, 2);
        //return std::pow(theta,2)/std::pow(std::exp(-theta*v/2)-std::exp(theta*v/2),2);
    };
    return u.unaryExpr(f);
}

// PDF
VecXd FrankBicop::pdf_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    MatXd t = u.unaryExpr([theta](const double v){ return std::exp(theta*v);});
    VecXd t1 = t.rowwise().prod();
    VecXd f = theta*(std::exp(theta)-1)*std::exp(theta)*t1;

    t1 = t1 - std::exp(theta)*(t.rowwise().sum() - VecXd::Ones(u.rows()));
    t1 = t1.array().square();

    f = f.cwiseQuotient(t1);
    return f;
}

// hfunction
VecXd FrankBicop::hfunc1_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    MatXd t = u.unaryExpr([theta](const double v){ return std::exp(theta*v);});
    VecXd t1 = t.rowwise().prod();
    VecXd f = std::exp(theta)*(t.col(1) - VecXd::Ones(u.rows()));

    t1 = - t1 + std::exp(theta)*(t.rowwise().sum() - VecXd::Ones(u.rows()));
    f = f.cwiseQuotient(t1);
    return f;
}

// inverse h-function
VecXd FrankBicop::hinv1_default(const MatXd& u)
{
    VecXd hinv = hinv1_num(u);
    return hinv;
}

VecXd FrankBicop::tau_to_parameters(const double& tau)
{
    VecXd tau2 = VecXd::Constant(1, std::fabs(tau));
    auto f = [&](const VecXd &v) {
        return VecXd::Constant(1, std::fabs(parameters_to_tau(v)));
    };
    return invert_f(tau2, f, -100+1e-6, 100);
}

double FrankBicop::parameters_to_tau(const VecXd& parameters)
{
    double par = parameters(0);
    double tau = 1 - 4/par;
    double d = debyen(std::fabs(par), 1) / std::fabs(par);
    if (par < 0)
        d = d - par/2;
    tau = tau + (4/par) * d;
    return tau;
}

VecXd FrankBicop::get_start_parameters(const double tau)
{
    return tau_to_parameters(tau);
}