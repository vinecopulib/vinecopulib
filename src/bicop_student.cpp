// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_student.hpp"
#include "tools_stats.hpp"
#include <cmath>
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

namespace vinecopulib
{
    StudentBicop::StudentBicop()
    {
        family_ = 2;
        family_name_ = "t";
        rotation_ = 0;
        association_direction_ = "both";
        parameters_ = VectorXd::Zero(2);
        parameters_(1) = 50.0;
        parameters_bounds_ = MatrixXd::Ones(2, 2);
        parameters_bounds_(0, 0) = -1.0;
        parameters_bounds_(1, 0) = 2.0;
        parameters_bounds_(1, 1) = 50.0;
    }

    StudentBicop::StudentBicop(const VectorXd& parameters)
    {
        StudentBicop();
        set_parameters(parameters);
    }

    StudentBicop::StudentBicop(const VectorXd& parameters, const int& rotation)
    {
        StudentBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

    VectorXd StudentBicop::pdf_default(const MatrixXd& u)
    {
        double rho = double(this->parameters_(0));
        double nu = double(this->parameters_(1));
        VectorXd f = VectorXd::Ones(u.rows());
        MatrixXd tmp = qt(u, nu);

        f = tmp.col(0).cwiseAbs2() + tmp.col(1).cwiseAbs2() - (2 * rho) * tmp.rowwise().prod();
        f /= nu * (1.0 - pow(rho, 2.0));
        f = f + VectorXd::Ones(u.rows());
        f = f.array().pow(-(nu + 2.0) / 2.0);
        f = f.cwiseQuotient(dt(tmp, nu).rowwise().prod());
        //f *= StableGammaDivision((nu + 2.0) / 2.0, nu / 2.0) / (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)));
        f *= boost::math::tgamma_ratio((nu + 2.0) / 2.0, nu / 2.0) / (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)));

        return f;
    }

    VectorXd StudentBicop::hfunc1_default(const MatrixXd& u)
    {
        double rho = double(this->parameters_(0));
        double nu = double(this->parameters_(1));
        VectorXd h = VectorXd::Ones(u.rows());
        MatrixXd tmp = qt(u, nu);
        h = nu * h + tmp.col(0).cwiseAbs2();
        h *= (1.0 - pow(rho, 2)) / (nu + 1.0);
        h = h.cwiseSqrt().cwiseInverse().cwiseProduct(tmp.col(1) - rho * tmp.col(0));
        h = pt(h, nu + 1.0);

        return h;
    }

    VectorXd StudentBicop::hinv1_default(const MatrixXd& u)
    {
        double rho = double(this->parameters_(0));
        double nu = double(this->parameters_(1));
        VectorXd hinv = VectorXd::Ones(u.rows());
        VectorXd tmp = u.col(1);
        VectorXd tmp2 = u.col(0);
        tmp = qt(tmp, nu + 1.0);
        tmp2 = qt(tmp2, nu);

        hinv = nu * hinv + tmp2.cwiseAbs2();
        hinv *= (1.0 - pow(rho, 2)) / (nu + 1.0);
        hinv = hinv.cwiseSqrt().cwiseProduct(tmp) + rho * tmp2;
        hinv = pt(hinv, nu );

        return hinv;
    }

    VectorXd StudentBicop::get_start_parameters(const double tau)
    {
        VectorXd parameters = tau_to_parameters(tau);
        parameters(1) = 5;
        return parameters;
    }
}
