//
// Created by Vatter Thibault on 13/12/16.
//

#ifndef __PARBICOP_H_INCLUDED__
#define __PARBICOP_H_INCLUDED__

#include "bicop.h"

class ParBiCop : public BiCop {

public:
    // constructors --------------------------------

    // default constructor (independence copula)
    ParBiCop();
    // construct with family and parameters
    ParBiCop(int family, double par, double par2);

    // getters and setters --------------------------------

    int getFamily();
    double getPar();
    double getPar2();
    double getNpars();
    void setFamily(int family);
    void setPar(double par);
    void setPar2(double par2);


    // calculate number of parameters
    int calculateNpars() ;

private:
    // copula specification
    int _family;     // copula family
    double _par;     // first copula parameter
    double _par2;    // second copula parameter
    int _npars;      // number of parameters
};


#endif //__PARBICOP_H_INCLUDED__
