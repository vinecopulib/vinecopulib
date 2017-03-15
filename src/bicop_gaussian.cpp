// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_gaussian.hpp"
#include "tools_stats.hpp"

namespace vinecopulib
{
    GaussianBicop::GaussianBicop()
    {
        family_ = 1;
        family_name_ = "Gaussian";
        rotation_ = 0;
        association_direction_ = "both";
        parameters_ = VectorXd::Zero(1);
        parameters_bounds_ = MatrixXd::Ones(1, 2);
        parameters_bounds_(0, 0) = -1;
    }

    GaussianBicop::GaussianBicop(const VectorXd& parameters)
    {
        GaussianBicop();
        set_parameters(parameters);
    }

    GaussianBicop::GaussianBicop(const VectorXd& parameters, const int& rotation)
    {
        GaussianBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

    VectorXd GaussianBicop::pdf_default(const MatrixXd& u)
    {
        // Inverse Cholesky of the correlation matrix
        double rho = double(this->parameters_(0));
        MatrixXd L = MatrixXd::Zero(2,2);
        L(0,0) = 1;
        L(1,1) = 1/sqrt(1.0-pow(rho,2.0));
        L(0,1) = -rho*L(1,1);

        // Compute copula density
        VectorXd f = VectorXd::Ones(u.rows());
        MatrixXd tmp = qnorm(u);
        f = f.cwiseQuotient(dnorm(tmp).rowwise().prod());
        tmp = tmp*L;
        f = f.cwiseProduct(dnorm(tmp).rowwise().prod());
        return f / sqrt(1.0-pow(rho,2.0));
    }

    VectorXd GaussianBicop::hfunc1_default(const MatrixXd& u)
    {
        double rho = double(this->parameters_(0));
        VectorXd h = VectorXd::Zero(u.rows());
        MatrixXd tmp = qnorm(u);
        h = (tmp.col(1) - rho * tmp.col(0)) / sqrt(1.0 - pow(rho, 2.0));
        return pnorm(h);
    }

    VectorXd GaussianBicop::hinv1_default(const MatrixXd& u)
    {
        double rho = double(this->parameters_(0));
        VectorXd hinv = VectorXd::Zero(u.rows());
        MatrixXd tmp = qnorm(u);
        hinv = tmp.col(1) * sqrt(1.0 - pow(rho, 2.0)) + rho * tmp.col(0);
        return pnorm(hinv);
    }

    VectorXd GaussianBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }
}
