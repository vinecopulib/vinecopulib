//
// Created by Vatter Thibault on 15/12/16.
//

#ifndef VINECOPLIB_JOE_BICOP_H
#define VINECOPLIB_JOE_BICOP_H


#include "archimedean_bicop.h"

class JoeBicop : public ArchimedeanBicop {

public:
    // constructor
    JoeBicop();
    JoeBicop(double theta);
    JoeBicop(double theta, int rotation);

    // generator, its inverse and derivatives for the archimedean copula
    VecXd generator(const VecXd &u);
    VecXd generator_inv(const VecXd &u);
    VecXd generator_derivative(const VecXd &u);
    VecXd generator_derivative2(const VecXd &u);

};


#endif //VINECOPLIB_JOE_BICOP_H
