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

#include "bicop_frank.hpp"


// constructor
FrankBicop::FrankBicop()
{
    family_ = 5;
    family_name_ = "Frank";
    rotation_ = 0;
    association_direction_ = "both";
    parameters_ = VecXd::Zero(1);
    parameters_bounds_ = MatXd::Zero(1, 2);
    parameters_bounds_(0, 0) = -200.0;
    parameters_bounds_(0, 1) = 200.0;
}

FrankBicop::FrankBicop(const VecXd& parameters)
{
    FrankBicop();
    set_parameters(parameters);
}

FrankBicop::FrankBicop(const VecXd& parameters, const int& rotation)
{
    FrankBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd FrankBicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-theta) * u;
    psi = psi.array().exp();
    psi = (psi - VecXd::Ones(psi.size())) / (exp(-theta) - 1);
    psi = (-1) * psi.array().log();
    return psi;
}
VecXd FrankBicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1) * u;
    psi = (exp(-theta) - 1) * psi.array().exp();
    psi = psi + VecXd::Ones(psi.size());
    psi = (-1 / theta) * psi.array().log();
    return psi;
}

VecXd FrankBicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = theta * u;
    psi = psi.array().exp();
    psi = VecXd::Ones(psi.size()) - psi;
    psi = theta * psi.cwiseInverse();
    return psi;
}

VecXd FrankBicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-(theta/2) * u).array().exp() - ((theta/2) * u).array().exp();
    psi = (theta*theta) * psi.cwiseInverse().array().square();
    return psi;
}

VecXd FrankBicop::hinv(const MatXd& u)
{
    VecXd hinv = hinv1_num(u);
    return hinv;
}

VecXd FrankBicop::tau_to_parameters(const double& tau)
{
    int br = 0, it = 0;
    double tol = 1e-12, xl = -100+1e-6, xh = 100, fl, fh, fm, par;

    fl = std::fabs(parameters_to_tau(VecXd::Constant(1, xl)) - tau);
    fh = std::fabs(parameters_to_tau(VecXd::Constant(1, xh)) - tau);
    if (fl <= tol) {
        par = xl;
        br = 1;
    }
    if (fh <= tol) {
        par = xh;
        br = 1;
    }

    while (!br) {
        par = (xh + xl) / 2.0;
        fm = parameters_to_tau(VecXd::Constant(1, par)) - tau;

        //stop if values become too close (avoid infinite loop)
        if (std::fabs(fm) <= tol) br = 1;
        if (std::fabs(xl - xh) <= tol) br = 1;

        if (tau < 0.0)
        {
            if (fm > 0.0) {
                xl = par;
                fh = fm;
            } else {
                xh = par;
                fl = fm;
            }
        } else
        {
            if (fm < 0.0) {
                xl = par;
                fh = fm;
            } else {
                xh = par;
                fl = fm;
            }
        }

        //stop if too many iterations are required (avoid infinite loop)
        ++it;
        if (it > 50)
            br = 1;
    }

    VecXd parameters = VecXd::Constant(1, par);
    return parameters;
}

double FrankBicop::parameters_to_tau(const VecXd& parameters)
{
    double par = parameters(0);
    double tau = 1 - 4/par;
    double d = gsl_sf_debye_1(std::fabs(par));
    if (par < 0)
        d = d - par/2;
    tau = tau + (4/par) * d;
    return tau;
}
