//
// Created by Vatter Thibault on 15/12/16.
//

#ifndef VINECOPLIB_FRANK_BICOP_H
#define VINECOPLIB_FRANK_BICOP_H

#include "archimedean_bicop.h"

class FrankBicop : public ArchimedeanBicop {

public:
    // constructor
    FrankBicop();
    FrankBicop(double theta);
    FrankBicop(double theta, int rotation);

    // generator, its inverse and derivatives for the archimedean copula
    VecXd generator(const VecXd &u);
    VecXd generator_inv(const VecXd &u);
    VecXd generator_derivative(const VecXd &u);
    VecXd generator_derivative2(const VecXd &u);

};


#endif //VINECOPLIB_FRANK_BICOP_H
