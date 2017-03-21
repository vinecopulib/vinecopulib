// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_bb6.hpp"
#include "tools_integration.hpp"

namespace vinecopulib
{
    Bb6Bicop::Bb6Bicop()
    {
        family_ = BicopFamily::bb6;
        rotation_ = 0;
        parameters_ = Eigen::VectorXd::Ones(2);
        parameters_bounds_ = Eigen::MatrixXd::Constant(2, 2, 200);
        parameters_bounds_(0, 0) = 1.0;
        parameters_bounds_(1, 0) = 1.0;
    }

    double Bb6Bicop::generator(const double& u)
    {
        return std::pow((-1)*std::log(1-std::pow(1-u,this->parameters_(0))), this->parameters_(1));
    }

    double Bb6Bicop::generator_inv(const double& u)
    {
        return 1-std::pow(1-std::exp((-1)*std::pow(u, 1/this->parameters_(1))), 1/this->parameters_(0));
    }

    double Bb6Bicop::generator_derivative(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        double res = delta * theta *std::pow((-1)*std::log(1-std::pow(1-u,theta)),delta-1);
        return res*std::pow(1-u,theta-1)/(std::pow(1-u,theta)-1);
    }

    double Bb6Bicop::generator_derivative2(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        double tmp = std::pow(1-u,theta);
        double res = std::pow((-1)*std::log(1-tmp),delta-2)*((delta-1)*theta*tmp-(tmp+theta-1)*std::log(1-tmp));
        return res*delta*theta*std::pow(1-u, theta-2)/std::pow(tmp - 1,2);
    }

    double Bb6Bicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        double theta = parameters(0);
        double delta = parameters(1);
        auto f = [&theta, &delta](const double& v) {
            return -4/(delta*theta)*std::log(1-std::pow(1-v,theta))*(1-v-std::pow(1-v,-theta)+std::pow(1-v,-theta)*v);
        };
        double tau = 1+tools_integration::integrate_zero_to_one(f);
        return flip_tau(tau);
    }

    Eigen::VectorXd Bb6Bicop::tau_to_parameters_default(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
