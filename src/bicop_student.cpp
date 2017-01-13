/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "bicop_student.hpp"

// constructor
StudentBicop::StudentBicop()
{
    family_ = 2;
    family_name_ = "t";
    rotation_ = 0;
    association_direction_ = "both";
    parameters_ = VecXd::Zero(2);
    parameters_(1) = 100.0;
    parameters_bounds_ = MatXd::Ones(2, 2);
    parameters_bounds_(0, 0) = -1.0;
    parameters_bounds_(1, 0) = 2.0;
    parameters_bounds_(1, 1) = 50.0;
}

StudentBicop::StudentBicop(const VecXd& parameters)
{
    StudentBicop();
    set_parameters(parameters);
}

StudentBicop::StudentBicop(const VecXd& parameters, const int& rotation)
{
    StudentBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

// PDF
VecXd StudentBicop::pdf_default(const MatXd& u)
{
    int j;
    double t1, t2;
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    VecXd f = VecXd::Zero(u.rows());

    for (j = 0; j < u.rows(); ++j) {
        t1 = gsl_cdf_tdist_Pinv(u(j, 0), nu);
        t2 = gsl_cdf_tdist_Pinv(u(j, 1), nu);
        f(j) = StableGammaDivision((nu + 2.0) / 2.0, nu / 2.0) /
        (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)) *
        gsl_ran_tdist_pdf(t1, nu) *
        gsl_ran_tdist_pdf(t2, nu)) *
        pow(1.0 + (pow(t1, 2.0) + pow(t2, 2.0) - 2.0 * rho * t1 * t2) /
        (nu * (1.0 - pow(rho, 2.0))),
        -(nu+2.0)/2.0
    );
}

return f;
}

// Student h-function
VecXd StudentBicop::hfunc1_default(const MatXd& u)
{
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    int j;
    double u1, u2, t1, t2, mu, sigma2;
    VecXd h = VecXd::Zero(u.rows());

    for(j = 0; j < u.rows(); ++j) {
        u1 = u(j, 1);
        u2 = u(j, 0);

        if((u1 == 0) | (u2 == 0)) {
            h(j) = 0; //else if (u2==1) h(j) = u1;
        } else {
            t1 = gsl_cdf_tdist_Pinv(u1,nu);
            t2 = gsl_cdf_tdist_Pinv(u2,nu);
            mu = rho * t2;
            sigma2 = ((nu + t2*t2) * (1.0 - rho*rho)) / (nu + 1.0);
            h(j) = gsl_cdf_tdist_P((t1 - mu) / sqrt(sigma2) , nu + 1.0);
        }
    }

    return h;
}


VecXd StudentBicop::hinv1_default(const MatXd& u)
{
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    int j;
    double u1, u2, t1, t2, mu, sigma2;
    VecXd hinv = VecXd::Zero(u.rows());

    for(j = 0; j < u.rows(); ++j) {
        u1 = u(j, 1);
        u2 = u(j, 0);
        t1 = gsl_cdf_tdist_Pinv(u1, nu + 1.0);
        t2 = gsl_cdf_tdist_Pinv(u2, nu);
        mu = rho * t2;
        sigma2 = (nu + t2*t2) * (1.0 - rho*rho) / (nu + 1.0);
        hinv(j) = gsl_cdf_tdist_P(sqrt(sigma2) * t1 + mu, nu);
    }

    return hinv;
}
