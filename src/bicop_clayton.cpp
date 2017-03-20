// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.


#include "bicop_clayton.hpp"
#include "tools_stl.hpp"

namespace vinecopulib
{
    ClaytonBicop::ClaytonBicop()
    {
        family_ = BicopFamily::clayton;
        rotation_ = 0;
        parameters_ = Eigen::VectorXd::Zero(1);
        parameters_bounds_ = Eigen::MatrixXd::Zero(1, 2);
        parameters_bounds_(0, 1) = 200.0;
    }

    ClaytonBicop::ClaytonBicop(const Eigen::VectorXd& parameters) :
        ClaytonBicop()
    {
        set_parameters(parameters);
    }

    ClaytonBicop::ClaytonBicop(const Eigen::VectorXd& parameters, const int& rotation) :
        ClaytonBicop(parameters)
    {
        set_rotation(rotation);
    }

    Eigen::VectorXd ClaytonBicop::generator(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [&theta](const double& v) {
            return (std::pow(v, -theta)-1)/theta;
        };
        return u.unaryExpr(f);
    }
    Eigen::VectorXd ClaytonBicop::generator_inv(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [&theta](const double& v) {
            return std::pow(1+theta*v, -1/theta);
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd ClaytonBicop::generator_derivative(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [&theta](const double& v) {
            return (-1)*std::pow(v, -1-theta);
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd ClaytonBicop::generator_derivative2(const Eigen::VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [&theta](const double& v) {
            return (1+theta)*std::pow(v, -2-theta);
        };
        return u.unaryExpr(f);
    }

    Eigen::VectorXd ClaytonBicop::hinv1_default(const Eigen::MatrixXd& u)
    {
        double theta = double(this->parameters_(0));
        Eigen::VectorXd hinv = u.col(0).array().pow(theta + 1.0);
        if (theta < 75) {
            hinv = u.col(1).cwiseProduct(hinv);
            hinv = hinv.array().pow(-theta/(theta + 1.0));
            Eigen::VectorXd x = u.col(0);
            x = x.array().pow(-theta);
            hinv = hinv - x + Eigen::VectorXd::Ones(x.size());
            hinv = hinv.array().pow(-1/theta);
        } else {
            hinv = hinv1_num(u);
        }
        return hinv;
    }

    Eigen::VectorXd ClaytonBicop::tau_to_parameters(const double& tau)
    {
        Eigen::VectorXd parameters(1);
        parameters(0) = 2 * std::fabs(tau) / (1 - std::fabs(tau));
        return parameters;
    }

    double ClaytonBicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        double tau =  parameters(0) / (2 + std::fabs(parameters(0)));
        return flip_tau(tau);
    }

    Eigen::VectorXd ClaytonBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }
}

/*// PDF
Eigen::VectorXd ClaytonBicop::pdf_default(const Eigen::MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    Eigen::VectorXd t1 = generator(u.col(0));
    Eigen::VectorXd t2 = generator(u.col(1));
    Eigen::VectorXd t = t1+t2;
    Eigen::VectorXd f = generator_inv(t);

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
Eigen::VectorXd ClaytonBicop::hfunc1_default(const Eigen::MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    Eigen::VectorXd t1 = generator(u.col(0));
    Eigen::VectorXd t2 = generator(u.col(1));
    Eigen::VectorXd t = t1+t2;
    Eigen::VectorXd f = generator_inv(t);
    f = f.array().pow(1+theta);

    t2 = u.col(0).array().pow(-1-theta);
    f = f.cwiseProduct(t2);
    return f;
}*/
