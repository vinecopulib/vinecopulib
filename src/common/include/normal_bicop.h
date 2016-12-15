//
// Created by Vatter Thibault on 14/12/16.
//

#ifndef VINECOPLIB_NORMAL_BICOP_H
#define VINECOPLIB_NORMAL_BICOP_H

#include "par_bicop.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
//#include <boost/math/distributions/normal.hpp>

class NormalBicop : public ParBicop {

public:
    // constructor
    NormalBicop();
    NormalBicop(double rho);

 /*   // hfunctions: the conditioning variable is put second
    VecXd hfunc1(const MatXd *u);
    VecXd hfunc2(const MatXd *u);

    // PDF
    VecXd pdf(const MatXd *u);*/
    // hfunctions: the conditioning variable is put second
    VecXd hfunc1(const MatXd &u);
    VecXd hfunc2(const MatXd &u);

    // PDF
    VecXd pdf(const MatXd &u);
};


#endif //VINECOPLIB_NORMAL_BICOP_H
