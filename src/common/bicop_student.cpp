//
// Created by Vatter Thibault on 21/12/16.
//

#include "include/bicop_student.hpp"

// constructor
StudentBicop::StudentBicop()
{
    family_ = 2;
    parameters_ = VecXd::Zero(2);
    parameters_(1) = 1e2;
}

StudentBicop::StudentBicop(double rho, double nu)
{
    family_ = 1;
    VecXd par = VecXd::Zero(1);
    par(0) = rho;
    par(1) = nu;
    parameters_ = par;
}

// Student h-function
VecXd StudentBicop::hfunc1(const MatXd &u)
{
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    int j;
    double u1, u2, t1, t2, mu, sigma2;
    VecXd h = VecXd::Zero(u.rows());

    for(j=0;j<u.rows();j++)
    {
        u1=u(j,1);
        u2=u(j,0);

        if((u1==0) | ( u2==0)) h(j) = 0;//else if (u2==1) h(j) = u1;
        else
        {
            t1 = gsl_cdf_tdist_Pinv(u1,nu);
            t2 = gsl_cdf_tdist_Pinv(u2,nu);
            mu = rho*t2; 
            sigma2 = ((nu+t2*t2)*(1.0-rho*(rho)))/(nu+1.0);
            h(j) = gsl_cdf_tdist_P((t1-mu)/sqrt(sigma2),nu+1.0);
        }
    }
    return(h);
}


VecXd StudentBicop::hinv1(const MatXd &u)
{
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    int j;
    double u1, u2, t1, t2, mu, sigma2;
    VecXd hinv = VecXd::Zero(u.rows());

    for(j=0;j<u.rows();j++)
    {
        u1 = u(j,1);
        u2 = u(j,0);
        t1 = gsl_cdf_tdist_Pinv(u1,nu+1.0);
        t2 = gsl_cdf_tdist_Pinv(u2,nu);
        mu = rho*t2;
        sigma2 = ((nu+(t2*t2))*(1.0-(rho*(rho))))/(nu+1.0);
        hinv(j) = gsl_cdf_tdist_P((sqrt(sigma2)*t1)+mu,nu);
    }
    return(hinv);

}

// PDF
VecXd StudentBicop::pdf(const MatXd &u)
{
    int j;
    double u1, u2, t1, t2;
    double rho = double(this->parameters_(0));
    double nu = double(this->parameters_(1));
    VecXd f = VecXd::Zero(u.rows());

    for(j=0;j<u.rows();j++)
    {
        u1=u(j,0);
        u2=u(j,1);
        t1 = gsl_cdf_tdist_Pinv(u1,nu); t2 = gsl_cdf_tdist_Pinv(u2,nu);
        f(j) = StableGammaDivision((nu+2.0)/2.0,nu/2.0)/(nu*M_PI*sqrt(1.0-pow(rho,2.0))*gsl_ran_tdist_pdf(t1,nu)*gsl_ran_tdist_pdf(t2,nu))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0);

    }
    return(f);
}

MatXd StudentBicop::get_bounds()
{
    MatXd bounds = MatXd::Ones(2,2);
    bounds(0,0) = -1;
    bounds(1,0) = 2;
    bounds(1,1) = 50;
    return(bounds);
}