// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_bb8.hpp"
#include "tools_integration.hpp"

namespace vinecopulib
{
    Bb8Bicop::Bb8Bicop()
    {
        family_ = BicopFamily::bb8;
        rotation_ = 0;
        parameters_ = Eigen::VectorXd::Ones(2);
        parameters_bounds_ = Eigen::MatrixXd::Constant(2, 2, 1);
        parameters_bounds_(0, 1) = 200.0;
        parameters_bounds_(1, 0) = 0.0;
    }

    Bb8Bicop::Bb8Bicop(const Eigen::VectorXd& parameters) : 
        Bb8Bicop()
    {
        set_parameters(parameters);
    }

    Bb8Bicop::Bb8Bicop(const Eigen::VectorXd& parameters, const int& rotation) :
        Bb8Bicop(parameters)
    {
        set_rotation(rotation);
    }

    Eigen::VectorXd Bb8Bicop::generator(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            return -std::log((1-std::pow(1-delta*v,theta))/(1-std::pow(1-delta,theta)));
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd Bb8Bicop::generator_inv(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            return (1-std::pow(1+std::exp(-v)*(std::pow(1-delta,theta)-1),1/theta))/delta;
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd Bb8Bicop::generator_derivative(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            return -delta*theta*std::pow(1-delta*v,theta-1)/(1-std::pow(1-delta*v,theta));
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd Bb8Bicop::generator_derivative2(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            double tmp = std::pow(1-delta*v,theta);
            return std::pow(delta,2)*theta*std::pow(1-delta*v,theta-2)*(theta-1+tmp)/std::pow(tmp-1,2);
        };
        return u.unaryExpr(f);
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
}
