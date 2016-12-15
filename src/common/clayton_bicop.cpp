//
// Created by Vatter Thibault on 15/12/16.
//

#include "include/clayton_bicop.h"

// constructor
ClaytonBicop::ClaytonBicop()
{
    family_ = 3;
    parameters_ = VecXd::Zero(1);
    rotation_ = 0;
}

ClaytonBicop::ClaytonBicop(double theta)
{
    family_ = 3;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = 0;
}

ClaytonBicop::ClaytonBicop(double theta, int rotation)
{
    family_ = 3;
    VecXd par = VecXd::Zero(1);
    par(0) = theta;
    parameters_ = par;
    rotation_ = rotation;
}

VecXd ClaytonBicop::generator(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = u.array().pow(-theta);
    psi = (psi-VecXd::Ones(psi.size()))/theta;
    return(psi);
}
VecXd ClaytonBicop::generator_inv(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (theta*u+VecXd::Ones(u.size()));
    psi = psi.array().pow(-1/theta);
    return(psi);
}

VecXd ClaytonBicop::generator_derivative(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (-1)*u.array().pow(-1-theta);
    return(psi);
}

VecXd ClaytonBicop::generator_derivative2(const VecXd &u)
{
    double theta = double(this->parameters_(0));
    VecXd psi = (1+theta)*u.array().pow(-2-theta);
    return(psi);
}