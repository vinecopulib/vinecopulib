//
// Created by Vatter Thibault on 15/12/16.
//

#include "include/frank_bicop.h"


// constructor
FrankBicop::FrankBicop()
{
    family_ = 5;
    parameters_ = VecXd::Zero(1);
    rotation_ = 0;
}

FrankBicop::FrankBicop(double theta)
{
    family_ = 5;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = 0;
}

FrankBicop::FrankBicop(double theta, int rotation)
{
    family_ = 5;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = rotation;
}

VecXd FrankBicop::generator(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-theta)*u;
    psi = psi.array().exp();
    psi = (psi - VecXd::Ones(psi.size()))/(exp(-theta)-1);
    psi = (-1)*psi.array().log();
    return(psi);
}
VecXd FrankBicop::generator_inv(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1)*u;
    psi = (exp(-theta)-1)*psi.array().exp();
    psi = psi+VecXd::Ones(psi.size());
    psi = (-1/theta)*psi.array().log();
    return(psi);
}

VecXd FrankBicop::generator_derivative(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = theta*u;
    psi = psi.array().exp();
    psi = VecXd::Ones(psi.size())-psi;
    psi = theta*psi.cwiseInverse();
    return(psi);
}

VecXd FrankBicop::generator_derivative2(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-(theta/2)*u).array().exp()-((theta/2)*u).array().exp();
    psi = (theta*theta)*psi.cwiseInverse().array().square();
    return(psi);
}