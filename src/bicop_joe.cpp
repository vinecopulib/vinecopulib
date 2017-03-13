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

#include "bicop_joe.hpp"

// constructor
JoeBicop::JoeBicop()
{
    family_ = 6;
    family_name_ = "Joe";
    rotation_ = 0;
    association_direction_ = "positive";
    parameters_ = VecXd::Ones(1);
    parameters_bounds_ = MatXd::Ones(1, 2);
    parameters_bounds_(0, 1) = 200.0;
}

JoeBicop::JoeBicop(const VecXd& parameters)
{
    JoeBicop();
    set_parameters(parameters);
}

JoeBicop::JoeBicop(const VecXd& parameters, const int& rotation)
{
    JoeBicop();
    set_parameters(parameters);
    set_rotation(rotation);
}

VecXd JoeBicop::generator(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (-1)*std::log(1-std::pow(1-v, theta));
    };
    return u.unaryExpr(f);
}

VecXd JoeBicop::generator_inv(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return 1-std::pow(1-std::exp(-v),1/theta);
    };
    return u.unaryExpr(f);
}

VecXd JoeBicop::generator_derivative(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return (-theta)*std::pow(1-v, theta-1)/(1-std::pow(1-v, theta));
    };
    return u.unaryExpr(f);
}

VecXd JoeBicop::generator_derivative2(const VecXd& u)
{
    double theta = double(this->parameters_(0));
    auto f = [theta](const double v) {
        return theta*(theta-1+std::pow(1-v, theta))*std::pow(1-v, theta-2)/std::pow(-1+std::pow(1-v, theta),2);
    };
    return u.unaryExpr(f);
}

// inverse h-function
VecXd JoeBicop::hinv1_default(const MatXd& u)
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
    VecXd tau2 = VecXd::Constant(1, std::fabs(tau));
    auto f = [&](const VecXd &v) {
        return VecXd::Constant(1, std::fabs(parameters_to_tau(v)));
    };
    return invert_f(tau2, f, 1+1e-6, 100);
}

double JoeBicop::parameters_to_tau(const VecXd& parameters)
{
    double par = parameters(0);
    double tau = 2 / par + 1;
    tau = boost::math::digamma(2.0) - boost::math::digamma(tau);
    tau = 1 + 2 * tau / (2 - par);
    if ((rotation_ == 90) | (rotation_ == 270))
        tau *= -1;
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
        if(std::isnan(pdf) || std::isnan(c21) ) { diff/=-2.; }  // added for de>=30
        else diff=(c21-*q)/pdf;
        v-=diff;
        int iter2 = 0;
        while((v<=0 || v>=1 || fabs(diff)>0.25) & (iter2 <20)) {++iter2; diff/=2.; v+=diff; }
    }
    return(v);
}

VecXd JoeBicop::get_start_parameters(const double tau)
{
    return tau_to_parameters(tau);
}

/*// PDF
VecXd JoeBicop::pdf_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));

    VecXd f = VecXd::Ones(u.rows());
    if (theta > 1+1e-6)
    {
        MatXd t = u.unaryExpr([theta](const double v){ return -1+std::pow(1-v,theta);});
        VecXd t1 = t.rowwise().prod();
        f = t1.unaryExpr([theta](const double v){ return theta-v;});

        t1 = t1.unaryExpr([theta](const double v){ return std::pow(1-v,-2+1/theta);});
        f = f.cwiseProduct(t1);

        t = u.unaryExpr([theta](const double v){ return std::pow(1-v,-1+theta);});
        t1 = t.rowwise().prod();
        f = f.cwiseProduct(t1);
    }
    return f;
}

// hfunction
VecXd JoeBicop::hfunc1_default(const MatXd& u)
{
    double theta = double(this->parameters_(0));
    MatXd t = u.unaryExpr([theta](const double v){ return -1+std::pow(1-v,theta);});
    VecXd t1 = t.rowwise().prod();
    VecXd f = t1.unaryExpr([theta](const double v){ return std::pow(1-v,-1+1/theta);});

    t1 = u.col(0).unaryExpr([theta](const double v){ return std::pow(1-v,-1+theta);});
    f = f.cwiseProduct(t1);

    t1 = u.col(1).unaryExpr([theta](const double v){ return 1-std::pow(1-v,theta);});
    f = f.cwiseProduct(t1);

    return f;
}*/
