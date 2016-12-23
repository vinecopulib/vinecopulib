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

#include "include/bicop_joe.hpp"

// constructor
JoeBicop::JoeBicop()
{
    family_ = 6;
    rotation_ = 0;
    parameters_ = VecXd::Ones(1);
    parameter_bounds_ = MatXd::Ones(1, 2);
    parameter_bounds_(0, 1) = 200.0;
}

JoeBicop::JoeBicop(const VecXd& parameters)
{
    family_ = 6;
    rotation_ = 0;
    parameters_ = parameters;
    parameter_bounds_ = MatXd::Ones(1, 2);
    parameter_bounds_(0, 1) = 200.0;
}

JoeBicop::JoeBicop(const VecXd& parameters, const int& rotation)
{
    family_ = 6;
    rotation_ = rotation;
    parameters_ = parameters;
    parameter_bounds_ = MatXd::Ones(1, 2);
    parameter_bounds_(0, 1) = 200.0;
}

VecXd JoeBicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = VecXd::Ones(u.size()) - u;
    psi = psi.array().pow(theta);
    psi = VecXd::Ones(psi.size()) - psi;
    psi = (-1) * psi.array().log();
    return psi;
}

VecXd JoeBicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1) * u;
    psi = psi.array().exp();
    psi = VecXd::Ones(psi.size()) - psi;
    psi = psi.array().pow(1/theta);
    psi = VecXd::Ones(psi.size()) - psi;
    return psi;
}

VecXd JoeBicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = VecXd::Ones(u.size()) - u;
    VecXd psi2 = psi.array().pow(theta);
    psi = psi.array().pow(-1+theta);
    psi2 = VecXd::Ones(psi.size()) - psi2;
    psi = (-theta) * psi.cwiseQuotient(psi2);
    return psi;
}

VecXd JoeBicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = VecXd::Ones(u.size()) - u;
    VecXd psi2 = psi.array().pow(theta);
    psi = psi.array().pow(-2+theta);
    psi2 = psi2 - VecXd::Ones(psi.size());
    psi = psi.cwiseQuotient(psi2.cwiseAbs2());
    psi = theta * psi.cwiseProduct(psi2 + theta * VecXd::Ones(psi.size()));
    return psi;
}

VecXd JoeBicop::hinv(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    double u1, u2;
    VecXd hinv = VecXd::Zero(u.rows());
    for (int j = 0; j < u.rows(); ++j) {
        u1 = u(j, 1);
        u2 = u(j, 0);
        hinv(j) = qcondjoe(&u1, &u2, &theta);
    }

    return hinv;
}

// link between Kendall's tau and the par_bicop parameter
VecXd JoeBicop::tau_to_parameters(const double& tau)
{
    int br = 0, it = 0;
    double tol = 1e-12, xl = -100 + 1e-6, xh = 100, fl, fh, fm, par;

    fl = fabs(parameters_to_tau(VecXd::Constant(1, xl)) - tau);
    fh = fabs(parameters_to_tau(VecXd::Constant(1, xh)) - tau);
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
        if (fabs(fm) <= tol) br = 1;
        if (fabs(xl - xh) <= tol) br = 1;

        if (fm < 0.0) {
            xl = par;
            fh = fm;
        } else {
            xh = par;
            fl = fm;
        }

        //stop if too many iterations are required (avoid infinite loop)
        ++it;
        if (it > 50)
            br = 1;
    }

    return VecXd::Constant(1, par);
}

double JoeBicop::parameters_to_tau(const VecXd& parameters)
{
    double par = parameters(0);
    double tau = 2 / par + 1;
    tau = gsl_sf_psi(2) - gsl_sf_psi(tau);
    tau = 1 + 2 * tau / (2 - par);
    return tau;
}

// This is copy&paste from the VineCopula package
double qcondjoe(double* q, double* u, double* de)
{
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t19,t23,t28,t31;
    double c21,pdf;
    int iter;
    double diff,v,de1,dtem,de1inv,tem;

    t1 = 1.0-*u;
    t2 = pow(t1,1.0*(*de));
    t7 = 1./(*de);
    t10 = t2*(*de);
    t11 = 1./t1;
    t19 = (*de)*(*de);
    de1=*de-1;  // may need better modification for large delta
    dtem=-de1/(1.+de1); de1inv=-1./de1;

    // v = 0.5 * (q+u); // starting guess

    // Use a better starting point based on reflected B4 copula
    // A good starting point is crucial when delta is large because
    //    C_{2|1} will be steep
    // C_{R,2|1}(v|u)=1-C_{2|1}(1-v|1-u),
    // C_{R,2|1}^{-1}(q|u)=1-C_{2|1}^{-1}(1-q|1-u)
    tem=pow(1.-*q,dtem)-1.;
    tem=tem*pow(1.-*u,-de1)+1.;
    v=pow(tem,de1inv); v=1.-v;
    diff=1; iter=0;
    while(fabs(diff)>1.e-6 && iter<20)
    { t3 = 1.-v;
        t4 = pow(t3,*de);
        t5 = t2*t4;
        t6 = t2+t4-t5;
        t8 = pow(t6,t7);
        t9 = t7*t8;
        t13 = t11*t4;
        t15 = -t10*t11+t10*t13;
        t16 = 1./t6;
        t23 = 1./t3;
        t28 = t6*t6;
        t31 = (-t4*(*de)*t23+t5*(*de)*t23)/t28*t15;
        c21 = -t9*t15*t16;
        pdf = -t8/t19*t31+t8*(*de)*t2*t13*t23*t16+t9*t31;
        iter++;
        if(isnan(pdf) || isnan(c21) ) { diff/=-2.; }  // added for de>=30
        else diff=(c21-*q)/pdf;
        v-=diff;
        while(v<=0 || v>=1 || fabs(diff)>0.25 ) { diff/=2.; v+=diff; }
    }
    return(v);
}
