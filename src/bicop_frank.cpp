// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_frank.hpp"
#include "tools_c.h"

namespace vinecopulib
{
    FrankBicop::FrankBicop()
    {
        family_ = 5;
        family_name_ = "Frank";
        rotation_ = 0;
        association_direction_ = "both";
        parameters_ = VectorXd::Zero(1);
        parameters_bounds_ = MatrixXd::Zero(1, 2);
        parameters_bounds_(0, 0) = -200.0;
        parameters_bounds_(0, 1) = 200.0;
    }

    FrankBicop::FrankBicop(const VectorXd& parameters)
    {
        FrankBicop();
        set_parameters(parameters);
    }

    FrankBicop::FrankBicop(const VectorXd& parameters, const int& rotation)
    {
        FrankBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

    VectorXd FrankBicop::generator(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return (-1)*std::log((std::exp(-theta*v)-1)/(std::exp(-theta)-1));
        };
        return u.unaryExpr(f);
    }
    VectorXd FrankBicop::generator_inv(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            //return (-1/theta)*std::log(1+std::exp(-theta*v)-std::exp(-v));
            return (-1/theta)*std::log(1+std::exp(-theta-v)-std::exp(-v));
        };
        return u.unaryExpr(f);
    }

    VectorXd FrankBicop::generator_derivative(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return theta/(1-std::exp(theta*v));
        };
        return u.unaryExpr(f);
    }

    VectorXd FrankBicop::generator_derivative2(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return std::pow(theta,2)/std::pow(std::exp(theta*v/2) - std::exp(-theta*v/2), 2);
        };
        return u.unaryExpr(f);
    }

    VectorXd FrankBicop::tau_to_parameters(const double& tau)
    {
        VectorXd tau2 = VectorXd::Constant(1, std::fabs(tau));
        auto f = [&](const VectorXd &v) {
            return VectorXd::Constant(1, std::fabs(parameters_to_tau(v)));
        };
        return invert_f(tau2, f, -100+1e-6, 100);
    }

    double FrankBicop::parameters_to_tau(const VectorXd& parameters)
    {
        double par = parameters(0);
        double tau = 1 - 4/par;
        double d = debyen(std::fabs(par), 1) / std::fabs(par);
        if (par < 0)
            d = d - par/2;
        tau = tau + (4/par) * d;
        return tau;
    }

    VectorXd FrankBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }
}
/*// PDF
VectorXd FrankBicop::pdf_default(const MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    MatrixXd t = u.unaryExpr([theta](const double v){ return std::exp(theta*v);});
    VectorXd t1 = t.rowwise().prod();
    VectorXd f = theta*(std::exp(theta)-1)*std::exp(theta)*t1;

    t1 = t1 - std::exp(theta)*(t.rowwise().sum() - VectorXd::Ones(u.rows()));
    t1 = t1.array().square();

    f = f.cwiseQuotient(t1);
    return f;
}

// hfunction
VectorXd FrankBicop::hfunc1_default(const MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    MatrixXd t = u.unaryExpr([theta](const double v){ return std::exp(theta*v);});
    VectorXd t1 = t.rowwise().prod();
    VectorXd f = std::exp(theta)*(t.col(1) - VectorXd::Ones(u.rows()));

    t1 = - t1 + std::exp(theta)*(t.rowwise().sum() - VectorXd::Ones(u.rows()));
    f = f.cwiseQuotient(t1);
    return f;
}*/
