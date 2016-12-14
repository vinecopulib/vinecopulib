//
// Created by Vatter Thibault on 13/12/16.
//

#ifndef __PARBICOP_H_INCLUDED__
#define __PARBICOP_H_INCLUDED__

#include "bicop.hpp"
#include <vector>

using namespace std;

class ParBiCop : public BiCop {

public:
    // constructors --------------------------------

    // default constructor (independence copula)
    ParBiCop();
    // construct with family and parameters
    ParBiCop(int new_family, vector<double> new_par);

    // getters and setters --------------------------------

    int get_family();
    vector<double> get_par();
    void set_family(int new_family);
    void set_par(vector<double> new_par);

    // number of parameters
    int calculate_npars();

private:
    // copula specification
    int family;     // copula family
    vector<double> par;     // first copula parameter
};


#endif //__PARBICOP_H_INCLUDED__
