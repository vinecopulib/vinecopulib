// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/bb1.hpp"
#include "misc/tools_integration.hpp"

namespace vinecopulib
{
    Bb1Bicop::Bb1Bicop()
    {
        family_ = BicopFamily::bb1;
        rotation_ = 0;
        parameters_ = Eigen::VectorXd(2);
        parameters_lower_bounds_ = Eigen::VectorXd(2);
        parameters_upper_bounds_ = Eigen::VectorXd(2);
        parameters_ << 0, 1;
        parameters_lower_bounds_ << 0, 1;
        parameters_upper_bounds_ << 200, 200;
    }


    double Bb1Bicop::generator(const double& u)
    {
        return std::pow(std::pow(u, -this->parameters_(0)) - 1, this->parameters_(1));
    }

    double Bb1Bicop::generator_inv(const double& u)
    {
        return std::pow(std::pow(u, 1/this->parameters_(1)) + 1, -1/this->parameters_(0));
    }

    double Bb1Bicop::generator_derivative(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        return -delta * theta * std::pow(u, -(1 + theta))*std::pow(std::pow(u, -theta) - 1, delta - 1);
    }

    double Bb1Bicop::generator_derivative2(const double& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        double res = delta * theta * std::pow(std::pow(u, -theta) - 1, delta) / std::pow(std::pow(u, theta) - 1, 2);
        return res * (1 + delta * theta - (1 + theta) * std::pow(u, theta)) / std::pow(u, 2);
    }

    double Bb1Bicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        double tau = 1-2/(parameters(1) * (parameters(0) + 2));
        return flip_tau(tau);
    }

    Eigen::MatrixXd Bb1Bicop::tau_to_parameters_default(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}


/*// PDF
double Bb1Bicop::pdf_default(const Eigen::MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));

    Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
    if (theta > 1e-6 || delta > 1+1e-6)
    {
        Eigen::VectorXd t1 = generator(u.col(0));
        Eigen::VectorXd t2 = generator(u.col(1));
        Eigen::VectorXd ones = f;
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
double Bb1Bicop::hfunc1_default(const Eigen::MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    Eigen::VectorXd t1 = generator(u.col(0));
    Eigen::VectorXd t2 = generator(u.col(1));
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(u.rows());
    Eigen::VectorXd f = ones;
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
