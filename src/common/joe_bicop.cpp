//
// Created by Vatter Thibault on 15/12/16.
//

#include "include/joe_bicop.h"

// constructor
JoeBicop::JoeBicop()
{
    family_ = 6;
    parameters_ = VecXd::Zero(1);
    rotation_ = 0;
}

JoeBicop::JoeBicop(double theta)
{
    family_ = 6;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = 0;
}

JoeBicop::JoeBicop(double theta, int rotation)
{
    family_ = 6;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = rotation;
}

VecXd JoeBicop::generator(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = VecXd::Ones(u.size())-u;
    psi = psi.array().pow(theta);
    psi = VecXd::Ones(psi.size())-psi;
    psi = (-1)*psi.array().log();
    return(psi);
}

VecXd JoeBicop::generator_inv(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1)*u;
    psi = psi.array().exp();
    psi = VecXd::Ones(psi.size())-psi;
    psi = psi.array().pow(1/theta);
    psi = VecXd::Ones(psi.size())-psi;
    return(psi);
}

VecXd JoeBicop::generator_derivative(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = VecXd::Ones(u.size())-u;
    VecXd psi2 = psi.array().pow(theta);
    psi = psi.array().pow(-1+theta);
    psi2 = VecXd::Ones(psi.size())-psi2;
    psi = (-theta)*psi.cwiseQuotient(psi2);
    return(psi);
}

VecXd JoeBicop::generator_derivative2(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = VecXd::Ones(u.size())-u;
    VecXd psi2 = psi.array().pow(theta);
    psi = psi.array().pow(-2+theta);
    psi2 = psi2-VecXd::Ones(psi.size());
    psi = psi.cwiseQuotient(psi2.cwiseAbs2());
    psi = theta*psi.cwiseProduct(psi2+theta*VecXd::Ones(psi.size()));
    return(psi);
}