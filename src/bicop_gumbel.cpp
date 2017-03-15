// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_gumbel.hpp"
#include <cmath>

namespace vinecopulib
{
    GumbelBicop::GumbelBicop()
    {
        family_ = 4;
        family_name_ = "Gumbel";
        rotation_ = 0;
        association_direction_ = "positive";
        parameters_ = VectorXd::Ones(1);
        parameters_bounds_ = MatrixXd::Ones(1, 2);
        parameters_bounds_(0, 1) = 200.0;
    }

    GumbelBicop::GumbelBicop(const VectorXd& parameters)
    {
        GumbelBicop();
        set_parameters(parameters);
    }

    GumbelBicop::GumbelBicop(const VectorXd& parameters, const int& rotation)
    {
        GumbelBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

    VectorXd GumbelBicop::generator(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return std::pow(std::log(1/v), theta);
        };
        return u.unaryExpr(f);
    }
    VectorXd GumbelBicop::generator_inv(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return std::exp(-std::pow(v,1/theta));
        };
        return u.unaryExpr(f);
    }

    VectorXd GumbelBicop::generator_derivative(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return std::pow(std::log(1/v),theta-1)*(-theta/v);
        };
        return u.unaryExpr(f);
    }

    VectorXd GumbelBicop::generator_derivative2(const VectorXd& u)
    {
        double theta = double(this->parameters_(0));
        auto f = [theta](const double v) {
            return (theta-1-std::log(v))*std::pow(std::log(1/v), theta-2)*(theta/std::pow(v,2));
        };
        return u.unaryExpr(f);
    }

    VectorXd GumbelBicop::hinv1_default(const MatrixXd& u)
    {
        double theta = double(this->parameters_(0));
        double u1, u2;
        VectorXd hinv = VectorXd::Zero(u.rows());
        for (int j = 0; j < u.rows(); ++j) {
            u1 = u(j, 1);
            u2 = u(j, 0);
            hinv(j) = qcondgum(&u1, &u2, &theta);
        }

        return hinv;
    }

    VectorXd GumbelBicop::tau_to_parameters(const double& tau)
    {
        return VectorXd::Constant(1, 1.0 / (1 - std::fabs(tau)));
    }

    double GumbelBicop::parameters_to_tau(const VectorXd& parameters)
    {
        double tau = (parameters(0) - 1) / parameters(0);
        if ((rotation_ == 90) | (rotation_ == 270))
            tau *= -1;

        return tau;
    }


    VectorXd GumbelBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }
}
// This is copy&paste from the VineCopula package
double qcondgum(double* q, double* u, double* de)
{
    double a,p,g,gp,z1,z2,con,de1,dif;
    double mxdif;
    int iter;

    p = 1-*q;
    z1 = -log(*u);
    con=log(1.-p)-z1+(1.-*de)*log(z1); de1=*de-1.;
    a=pow(2.*pow(z1,*de),1./(*de));
    mxdif=1; iter=0;
    dif=.1;  // needed in case first step leads to NaN
    while ((mxdif > 1.e-6) && (iter < 20))
    {
        g=a+de1*log(a)+con;
        gp=1.+de1/a;
        if (std::isnan(g) || std::isnan(gp) || std::isnan(g/gp) ) {
            // added for de>50
            dif/=-2.;
        } else {
            dif=g/gp;
        }
        a-=dif; iter++;
        int it = 0;
        while ((a <= z1) && (it < 20)) {
            dif /= 2.;
            a += dif;
            ++it;
        }
        mxdif=fabs(dif);
    }
    z2=pow(pow(a,*de)-pow(z1,*de),1./(*de));
    return(exp(-z2));
}
/*// PDF
VectorXd GumbelBicop::pdf_default(const MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    VectorXd t1 = generator(u.col(0));
    VectorXd t2 = generator(u.col(1));
    VectorXd t = t1+t2;
    t1 = t1.array().pow((theta-1)/theta);
    t2 = t2.array().pow((theta-1)/theta);
    VectorXd f = generator_inv(t);
    t = t.unaryExpr([theta](const double v){ return (-1+theta+std::pow(v, 1/theta)) * std::pow(v,-2+1/theta);});
    f = f.cwiseProduct(t).cwiseProduct(t1).cwiseProduct(t2).cwiseQuotient(u.rowwise().prod());
    return f;
}

// hfunction
VectorXd GumbelBicop::hfunc1_default(const MatrixXd& u)
{
    double theta = double(this->parameters_(0));
    VectorXd t1 = generator(u.col(1));
    VectorXd t2 = generator(u.col(0));
    VectorXd t = t1+t2;
    t2 = t2.array().pow((theta-1)/theta);
    VectorXd f = generator_inv(t);
    t = t.array().pow(1/theta-1);
    f = f.cwiseProduct(t).cwiseProduct(t2).cwiseQuotient(u.col(0));
    return f;
}*/
