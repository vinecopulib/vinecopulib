// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/bb8.hpp"
#include "misc/tools_integration.hpp"

namespace vinecopulib
{
    Bb8Bicop::Bb8Bicop()
    {
        family_ = BicopFamily::bb8;
        rotation_ = 0;
        parameters_ = Eigen::VectorXd(2);
        parameters_lower_bounds_ = Eigen::VectorXd(2);
        parameters_upper_bounds_ = Eigen::VectorXd(2);
        parameters_ << 1, 1;
        parameters_lower_bounds_ << 1, 0;
        parameters_upper_bounds_ << 200, 1;
    }

    double Bb8Bicop::generator(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        return -std::log((1-std::pow(1-delta*u,theta))/(1-std::pow(1-delta,theta)));
    }

    double Bb8Bicop::generator_inv(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        return (1-std::pow(1+std::exp(-u)*(std::pow(1-delta,theta)-1),1/theta))/delta;
    }

    double Bb8Bicop::generator_derivative(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        return -delta*theta*std::pow(1-delta*u,theta-1)/(1-std::pow(1-delta*u,theta));
    }

    double Bb8Bicop::generator_derivative2(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        double tmp = std::pow(1-delta*u,theta);
        return std::pow(delta,2)*theta*std::pow(1-delta*u,theta-2)*(theta-1+tmp)/std::pow(tmp-1,2);
    }

    double Bb8Bicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        double theta = parameters(0);
        double delta = parameters(1);
        auto f = [theta, delta](const double t) {
            double tmp = std::pow(1-t*delta,theta);
            return std::log((tmp-1)/(std::pow(1-delta,theta)-1))*(1-t*delta-std::pow(1-t*delta,1-theta));
        };
        double tau = 1-4/(delta*theta)*tools_integration::integrate_zero_to_one(f);
        return flip_tau(tau);
    }

    Eigen::MatrixXd Bb8Bicop::tau_to_parameters_default(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
