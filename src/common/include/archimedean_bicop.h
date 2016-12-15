//
// Created by Vatter Thibault on 15/12/16.
//

#ifndef VINECOPLIB_ARCHIMEDEAN_BICOP_H
#define VINECOPLIB_ARCHIMEDEAN_BICOP_H

#include "par_bicop.h"

class ArchimedeanBicop : public ParBicop {

public:
    // hfunctions: the conditioning variable is put second
    VecXd hfunc1(const MatXd &u);
    VecXd hfunc2(const MatXd &u);

    // PDF
    VecXd pdf(const MatXd &u);

    // generator, its inverse and derivatives for the archimedean copula
    virtual VecXd generator(const VecXd &u) = 0;
    virtual VecXd generator_inv(const VecXd &u) = 0;
    virtual VecXd generator_derivative(const VecXd &u) = 0;
    virtual VecXd generator_derivative2(const VecXd &u) = 0;

protected:
    int rotation_;
};


#endif //VINECOPLIB_ARCHIMEDEAN_BICOP_H
