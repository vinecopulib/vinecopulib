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
        family_ = BicopFamily::frank;
        rotation_ = 0;
        parameters_ = {Eigen::VectorXd(1)};
        parameters_lower_bounds_ = {Eigen::VectorXd(1)};
        parameters_upper_bounds_ = {Eigen::VectorXd(1)};
        parameters_[0] << 0;
        parameters_lower_bounds_[0] << -200;
        parameters_upper_bounds_[0] << 200;
    }

    Eigen::VectorXd FrankBicop::generator(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_[0](0));
        auto f = [theta](const double v) {
            return (-1)*std::log((std::exp(-theta*v)-1)/(std::exp(-theta)-1));
        };
        return u.unaryExpr(f);
    }
    Eigen::VectorXd FrankBicop::generator_inv(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_[0](0));
        auto f = [theta](const double v) {
            //return (-1/theta)*std::log(1+std::exp(-theta*v)-std::exp(-v));
            return (-1/theta)*std::log(1+std::exp(-theta-v)-std::exp(-v));
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd FrankBicop::generator_derivative(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_[0](0));
        auto f = [theta](const double v) {
            return theta/(1-std::exp(theta*v));
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd FrankBicop::generator_derivative2(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_[0](0));
        auto f = [theta](const double v) {
            return std::pow(theta,2)/std::pow(std::exp(theta*v/2) - std::exp(-theta*v/2), 2);
        };
        return u.unaryExpr(f);
    }

    std::vector<Eigen::MatrixXd> FrankBicop::tau_to_parameters(const double& tau)
    {
        Eigen::VectorXd tau2 = Eigen::VectorXd::Constant(1, std::fabs(tau));
        auto f = [&](const Eigen::VectorXd &v) {
            return Eigen::VectorXd::Constant(1, std::fabs(parameters_to_tau({v})));
        };
        return {invert_f(tau2, f, -100 + 1e-6, 100)};
    }

    double FrankBicop::parameters_to_tau(const std::vector<Eigen::MatrixXd>& parameters)
    {
        double par = parameters[0](0);
        double tau = 1 - 4 / par;
        double deb = debyen(std::fabs(par), 1) / std::fabs(par);
        if (par < 0) {
            deb = deb - par / 2;
        }
        return tau + (4 / par) * deb;
    }

    std::vector<Eigen::MatrixXd> FrankBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }
}
/*// PDF
Eigen::VectorXd FrankBicop::pdf_default(const Eigen::MatrixXd& u)
{
    double theta = double(this->parameters_[0](0));
    Eigen::MatrixXd t = u.unaryExpr([theta](const double v){ return std::exp(theta*v);});
    Eigen::VectorXd t1 = t.rowwise().prod();
    Eigen::VectorXd f = theta*(std::exp(theta)-1)*std::exp(theta)*t1;

    t1 = t1 - std::exp(theta)*(t.rowwise().sum() - Eigen::VectorXd::Ones(u.rows()));
    t1 = t1.array().square();

    f = f.cwiseQuotient(t1);
    return f;
}

// hfunction
Eigen::VectorXd FrankBicop::hfunc1_default(const Eigen::MatrixXd& u)
{
    double theta = double(this->parameters_[0](0));
    Eigen::MatrixXd t = u.unaryExpr([theta](const double v){ return std::exp(theta*v);});
    Eigen::VectorXd t1 = t.rowwise().prod();
    Eigen::VectorXd f = std::exp(theta)*(t.col(1) - Eigen::VectorXd::Ones(u.rows()));

    t1 = - t1 + std::exp(theta)*(t.rowwise().sum() - Eigen::VectorXd::Ones(u.rows()));
    f = f.cwiseQuotient(t1);
    return f;
}*/
