//
// Created by Vatter Thibault on 15/12/16.
//

#include "include/gumbel_bicop.h"

// constructor
GumbelBicop::GumbelBicop()
{
    family_ = 4;
    parameters_ = VecXd::Zero(1);
    rotation_ = 0;
}

GumbelBicop::GumbelBicop(double theta)
{
    family_ = 4;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = 0;
}

GumbelBicop::GumbelBicop(double theta, int rotation)
{
    family_ = 4;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = rotation;
}

VecXd GumbelBicop::generator(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = u;
    psi = psi.cwiseInverse().array().log().pow(theta);
    return(psi);
}
VecXd GumbelBicop::generator_inv(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = u;
    psi = (-1.0)*psi.array().pow(1/theta);
    psi = psi.array().exp();
    return(psi);
}

VecXd GumbelBicop::generator_derivative(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = u;
    psi = psi.cwiseInverse().array().log().pow(theta-1.0);
    psi = (-theta)*psi.cwiseQuotient(u);
    return(psi);
}

VecXd GumbelBicop::generator_derivative2(const VecXd &u)
{
    double theta = double(this->parameters_(0));;
    VecXd psi_ilog = u.cwiseInverse().array().log();
    VecXd psi = (theta-1)*VecXd::Ones(u.size())+psi_ilog;
    psi_ilog = psi_ilog.array().pow(theta-2.0);
    psi = theta*psi.cwiseProduct(psi_ilog);
    psi_ilog = u.array().square();
    psi = psi.cwiseQuotient(psi_ilog);
    return(psi);
}