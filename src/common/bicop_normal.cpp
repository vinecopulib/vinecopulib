/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecoplib.

    vinecoplib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecoplib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "include/bicop_normal.hpp"

// constructor
NormalBicop::NormalBicop()
{
    family_ = 1;
    parameters_ = VecXd::Zero(1);
}

NormalBicop::NormalBicop(double rho)
{
    family_ = 1;
    VecXd par = VecXd::Zero(1);
    par(0) = rho;
    parameters_ = par;
}

// Normal h-function
VecXd NormalBicop::hfunc1(const MatXd &u)
{
    double rho = double(this->parameters_(0));
    int j;
    double u1, u2, t1, t2;
    VecXd h = VecXd::Zero(u.rows());
    //boost::math::normal dist(0,1);

    for(j=0;j<u.rows();j++)
    {
        u1=u(j,0);
        u2=u(j,1);

        if((u1==0) | ( u2==0)) h(j) = 0;//else if (u2==1) h(j) = u1;
        else
        {
            t1 = gsl_cdf_ugaussian_Pinv(u1); t2 = gsl_cdf_ugaussian_Pinv(u2);
            //t1 = boost::math::quantile(dist, u1);
            //t2 = boost::math::quantile(dist, u2);
            h(j) = (t2 - rho*t1)/sqrt(1.0-pow(rho,2.0));
            if (std::isfinite(h(j)))
                h(j) = gsl_cdf_ugaussian_P(h(j));
                //h(j) = boost::math::cdf(dist, h(j));
            else if ((t2 - rho*t1) < 0)
                h(j) = 0;
            else
                h(j) = 1;
        }
    }
    return(h);
}


VecXd NormalBicop::hinv1(const MatXd &u)
{
    double rho = double(this->parameters_(0));
    VecXd hinv = VecXd::Zero(u.rows());

    for(int j=0;j<u.rows();j++)
    {
        hinv(j) = gsl_cdf_ugaussian_P(gsl_cdf_ugaussian_Pinv(u(j,1))*sqrt(1.0-pow(rho,2.0))+rho*gsl_cdf_ugaussian_Pinv(u(j,0)));
    }
    return(hinv);

}

// PDF
VecXd NormalBicop::pdf(const MatXd &u)
{
    int j;
    double u1, u2, t1, t2;
    double rho = double(this->parameters_(0));
    VecXd f = VecXd::Zero(u.rows());
    //boost::math::normal dist(0,1);

    for(j=0;j<u.rows();j++)
    {
        u1=u(j,0);
        u2=u(j,1);
        t1 = gsl_cdf_ugaussian_Pinv(u1); t2 = gsl_cdf_ugaussian_Pinv(u2);
        f(j) = gsl_ran_bivariate_gaussian_pdf(t1, t2, 1.0, 1.0, rho)/(gsl_ran_ugaussian_pdf(t1)*gsl_ran_ugaussian_pdf(t2));
        //t1 = boost::math::quantile(dist, u1);
        //t2 = boost::math::quantile(dist, u2);
        //f(j) = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
    }
    return(f);
}

MatXd NormalBicop::get_bounds()
{
    MatXd bounds = MatXd::Ones(1,2);
    bounds(0,0) = -1;
    return(bounds);
}