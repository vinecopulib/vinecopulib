// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_bb1.hpp"
#include "tools_integration.hpp"

namespace vinecopulib
{
    Bb1Bicop::Bb1Bicop()
    {
        family_ = BicopFamily::BB1;
        family_name_ = "BB1";
        rotation_ = 0;
        parameters_ = Eigen::VectorXd::Zero(2);
        parameters_(1) = 1;
        parameters_bounds_ = Eigen::MatrixXd::Constant(2, 2, 200);
        parameters_bounds_(0, 0) = 0.0;
        parameters_bounds_(1, 0) = 1.0;
    }

    Bb1Bicop::Bb1Bicop(const Eigen::VectorXd& parameters) : 
        Bb1Bicop()
    {
        set_parameters(parameters);
    }

    Bb1Bicop::Bb1Bicop(const Eigen::VectorXd& parameters, const int& rotation) :
        Bb1Bicop(parameters)
    {
        set_rotation(rotation);
    }

    Eigen::VectorXd Bb1Bicop::generator(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            return std::pow(std::pow(v, -theta) - 1, delta);
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd Bb1Bicop::generator_inv(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            return std::pow(std::pow(v, 1/delta) + 1, -1/theta);
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd Bb1Bicop::generator_derivative(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            return -delta * theta * std::pow(v, -(1 + theta))*std::pow(std::pow(v, -theta) - 1, delta - 1);
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd Bb1Bicop::generator_derivative2(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        double delta = double(this->parameters_(1));
        auto f = [theta, delta](const double v) {
            double res = delta * theta * std::pow(std::pow(v, -theta) - 1, delta) / std::pow(std::pow(v, theta) - 1, 2);
            return res * (1 + delta * theta - (1 + theta) * std::pow(v, theta)) / std::pow(v, 2);
        };
        return u.unaryExpr(f);
    }

    double Bb1Bicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        double tau = 1-2/(parameters(1) * (parameters(0) + 2));
        if ((rotation_ == 90) | (rotation_ == 270))
            tau *= -1;
        return tau;
    }
}

/*// PDF
Eigen::VectorXd Bb1Bicop::pdf_default(const Eigen::MatrixXd& u)
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
Eigen::VectorXd Bb1Bicop::hfunc1_default(const Eigen::MatrixXd& u)
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
