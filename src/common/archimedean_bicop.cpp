//
// Created by Vatter Thibault on 15/12/16.
//

#include "include/archimedean_bicop.h"

VecXd ArchimedeanBicop::hfunc1(const MatXd &u)
{
    VecXd v = generator(u.col(0))+generator(u.col(1));
    VecXd h = generator_derivative(u.col(0)).cwiseQuotient(generator_derivative(generator_inv(v)));
    return(h);
}

VecXd ArchimedeanBicop::hfunc2(const MatXd &u)
{
    MatXd v = u;
    v.col(0).swap(v.col(1));
    VecXd h = hfunc1(v);
    return(h);
}

// PDF
VecXd ArchimedeanBicop::pdf(const MatXd &u)
{
    VecXd v = generator(u.col(0))+generator(u.col(1));
    VecXd f = generator_derivative(u.col(0)).cwiseProduct(generator_derivative(u.col(1)));
    VecXd numerator = generator_derivative2(generator_inv(v));
    VecXd denominator = generator_derivative(generator_inv(v)).array().pow(3.0);
    f = (-1)*f.cwiseProduct(numerator);
    f = f.cwiseQuotient(denominator);
    return(f);
}
