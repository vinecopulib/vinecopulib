//
// Created by Vatter Thibault on 13/12/16.
//

#ifndef VINECOPLIB_PAR_BICOP_H_INCLUDED__
#define VINECOPLIB_PAR_BICOP_H_INCLUDED__

#include "bicop.h"

class ParBicop : public Bicop {

public:

    // getters and setters --------------------------------
    int get_family();
    VecXd get_parameters();
    void set_family(int family);
    void set_parameters(VecXd parameters);

    // other methods --------------------------------
    // number of parameters
    int calculate_npars();

protected:
    // copula specification
    int family_;                // copula family
    VecXd parameters_;       // first copula parameter
};


#endif //VINECOPLIB_PAR_BICOP_H_INCLUDED__
