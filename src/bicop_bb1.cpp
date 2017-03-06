//
// Created by Vatter Thibault on 22/01/17.
//

#include "bicop_bb1.hpp"


// constructor
Bb1Bicop::Bb1Bicop()
{
    family_ = 7;
    family_name_ = "Bb1";
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Zero(2);
    parameters_(1) = 1;
    parameters_bounds_ = MatXd::Constant(2, 2, 200);
    parameters_bounds_(0, 0) = 0.0;
    parameters_bounds_(1, 0) = 1.0;
}

Bb1Bicop::Bb1Bicop(const VecXd& parameters)
{
    Bb1Bicop();
    set_parameters(parameters);
}

Bb1Bicop::Bb1Bicop(const VecXd& parameters, const int& rotation)
{
    Bb1Bicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd Bb1Bicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return std::pow(std::pow(v, -theta) - 1, delta);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return std::pow(std::pow(v, 1/delta) + 1, -1/theta);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        return -delta * theta * std::pow(v, -(1 + theta))*std::pow(std::pow(v, -theta) - 1, delta - 1);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    double delta = double(this->parameters_(1));
    auto f = [theta, delta](const double v) {
        double res = delta * theta * std::pow(std::pow(v, -theta) - 1, delta) / std::pow(std::pow(v, theta) - 1, 2);
        return res * (1 + delta * theta - (1 + theta) * std::pow(v, theta)) / std::pow(v, 2);
    };
    return u.unaryExpr(f);
}

VecXd Bb1Bicop::hinv(const MatXd& u)
{
    VecXd hinv = hinv1_num(u);
    return hinv;
}

// link between Kendall's tau and the par_bicop parameter
VecXd Bb1Bicop::tau_to_parameters(const double& tau)
{
    double delta = double(this->parameters_(1));
    VecXd tau2 = VecXd::Constant(1, std::fabs(tau));
    auto f = [this, delta](const VecXd &v) {
        VecXd par = VecXd::Constant(2, delta);
        par(0) = v(0);
        return VecXd::Constant(1, std::fabs(this->parameters_to_tau(par)));
    };
    VecXd par = VecXd::Constant(2, delta);
    par(0) = invert_f(tau2, f, 1+1e-6, 100)(0);
    return par;
}

double Bb1Bicop::parameters_to_tau(const VecXd& parameters)
{
    double tau = 1-2/(parameters(1) * (parameters(0) + 2));
    if ((rotation_ == 90) | (rotation_ == 270))
        tau *= -1;
    return tau;
}

VecXd Bb1Bicop::get_start_parameters(const double tau)
{
    MatXd bounds = this->get_parameters_bounds();
    VecXd parameters = bounds.col(0) + VecXd::Constant(2, 0.1);
    return parameters;
}