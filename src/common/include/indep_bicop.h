//
// Created by Vatter Thibault on 14/12/16.
//

#ifndef VINECOPLIB_INDEP_BICOP_H
#define VINECOPLIB_INDEP_BICOP_H

#include "par_bicop.h"

class IndepBicop : public ParBicop {

public:
    // constructor
    IndepBicop();

    // hfunctions: the conditioning variable is put second
    VecXd hfunc1(const MatXd &u);
    VecXd hfunc2(const MatXd &u);

    // PDF
    VecXd pdf(const MatXd &u);
};


#endif //VINECOPLIB_INDEP_BICOP_H
