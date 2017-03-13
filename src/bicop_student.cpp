/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
    parameters_(1) = 50.0;
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
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    VecXd f = VecXd::Ones(u.rows());
    MatXd tmp = qt(u, nu);

    f = tmp.col(0).cwiseAbs2() + tmp.col(1).cwiseAbs2() - (2 * rho) * tmp.rowwise().prod();
    f /= nu * (1.0 - pow(rho, 2.0));
    f = f + VecXd::Ones(u.rows());
    f = f.array().pow(-(nu + 2.0) / 2.0);
    f = f.cwiseQuotient(dt(tmp, nu).rowwise().prod());
    //f *= StableGammaDivision((nu + 2.0) / 2.0, nu / 2.0) / (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)));
    f *= boost::math::tgamma_ratio((nu + 2.0) / 2.0, nu / 2.0) / (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)));

    return f;
}

// Student h-function
VecXd StudentBicop::hfunc1_default(const MatXd& u)
{
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    VecXd h = VecXd::Ones(u.rows());
    MatXd tmp = qt(u, nu);
    h = nu * h + tmp.col(0).cwiseAbs2();
    h *= (1.0 - pow(rho, 2)) / (nu + 1.0);
    h = h.cwiseSqrt().cwiseInverse().cwiseProduct(tmp.col(1) - rho * tmp.col(0));
    h = pt(h, nu + 1.0);

    return h;
}


VecXd StudentBicop::hinv1_default(const MatXd& u)
{
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    VecXd hinv = VecXd::Ones(u.rows());
    VecXd tmp = u.col(1);
    VecXd tmp2 = u.col(0);
    tmp = qt(tmp, nu + 1.0);
    tmp2 = qt(tmp2, nu);

    hinv = nu * hinv + tmp2.cwiseAbs2();
    hinv *= (1.0 - pow(rho, 2)) / (nu + 1.0);
    hinv = hinv.cwiseSqrt().cwiseProduct(tmp) + rho * tmp2;
    hinv = pt(hinv, nu );

    return hinv;
}

VecXd StudentBicop::get_start_parameters(const double tau)
{
    VecXd parameters = tau_to_parameters(tau);
    parameters(1) = 5;
    return parameters;
}
