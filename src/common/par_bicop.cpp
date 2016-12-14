//
// Created by Vatter Thibault on 13/12/16.
//

#include "include/par_bicop.hpp"

namespace {
    typedef Eigen::VectorXd dvec;
}


// constructors --------------------------------

// default constructor (independence copula)
ParBiCop::ParBiCop() :
        family(0),
        par(0)
{}

// construct with family and parameters
ParBiCop::ParBiCop(int new_family, vector<double> new_par)
{
    family = new_family;
    par = new_par;
}

// getters --------------------------------

int ParBiCop::get_family()
{
    return this->family;
}

vector<double> ParBiCop::get_par()
{
    return this->par;
}

// setters --------------------------------

void ParBiCop::set_family(int new_family)
{
    family = new_family;
}

void ParBiCop::set_par(vector<double> new_par)
{
    par = new_par;
}

// Family-specific functions --------------------------------

// calculate number of parameters
int ParBiCop::calculate_npars()
{
    // indepence copula has no parameters
    if (family == 0)
        return(0);

    // otherwise, return length of parameter vector
    return(int(par.size()));
}

/*
// hfunctions: the conditioning variable is put second
dvec ParBiCop::hfunc1(dvec &u1, dvec &u2)
{
    int n = u1.size();
    dvec out(n);

    //Hfunc1(&family, &n, &u2[0], &u1[0], &par[0], &out[0]);

    return out;
}

dvec ParBiCop::hfunc2(dvec &u1, dvec &u2)
{
    int n = u1.size();
    dvec out(n);

    //Hfunc2(&family, &n, &u1[0], &u2[0], &par[0], &out[0]);

    return out;
}

// PDF
dvec ParBiCop::pdf(dvec &u1, dvec &u2)
{
    int n = u1.size();
    dvec out(n);

    //PDF_seperate(&family, &n, &u1[0], &u2[0], &par, &out[0]);

    return out;
};

// fit statistics
double ParBiCop::loglik(dvec &u1, dvec &u2)
{
    int n = u1.size();
    double out;

    LL_mod2(&family, &n, &u1[0], &u2[0], &par, &out);

    return out;
}

double ParBiCop::aic(dvec &u1, dvec &u2)
{
    int n = u1.size();
    double out;

    //LL_mod2(&family, &n, &u1[0], &u2[0], &par, &out);
    out = -2 * out + 2 * calculate_npars();

    return out;
}

double ParBiCop::bic(dvec &u1, dvec &u2)
{
    int n = u1.size();
    double out;

    //LL_mod2(&family, &n, &u1[0], &u2[0], &par, &out);
    out = -2 * out + log(n) * calculate_npars();

    return out;
}*/
