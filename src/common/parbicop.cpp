//
// Created by Vatter Thibault on 13/12/16.
//

#include "include/parbicop.hpp"
#include <set>
#include <boost/assign/list_of.hpp>
#include <eigen3/Eigen/Dense>

namespace {
    typedef Eigen::VectorXd dvec;
}


// constructors --------------------------------

// default constructor (independence copula)
ParBiCop::ParBiCop() :
        _family(0),
        _par(0),
        _par2(0),
        _npars(0)
{}

// construct with family and parameters
ParBiCop::ParBiCop(int family, double par, double par2)
{
    _family = family;
    _par = par;
    _par2 = par2;
    _npars = calculateNpars();
}

// getters --------------------------------

int ParBiCop::getFamily()
{
    return this->_family;
}

double ParBiCop::getPar()
{
    return this->_par;
}

double ParBiCop::getPar2()
{
    return this->_par2;
}

double ParBiCop::getNpars()
{
    return this->_npars;
}

// setters --------------------------------

void ParBiCop::setFamily(int family)
{
    _family = family;
    _npars = calculateNpars();
}

void ParBiCop::setPar(double par)
{
    _par = par;
}

void ParBiCop::setPar2(double par2)
{
    _par2 = par2;
}

// Family-specific functions --------------------------------

// calculate number of parameters
int ParBiCop::calculateNpars()
{
    // indepence copula has no parameters
    if (_family == 0)
        return(0);

    // check if it is a one-parameter family
    std::set<int> values = boost::assign::list_of(1)(3)(4)(5)(6)(13)(14)(16)
            (23)(24)(26)(33)(34)(36);
    if (values.count(_family))
        return(1);

    // if not, it must be a two-parameter family
    return(2);
}

// hfunctions: the conditioning variable is put second
dvec ParBiCop::hFunc1(dvec &u1, dvec &u2)
{
    int n = u1.size();
    dvec out(n);

    Hfunc1(&_family, &n, &u2[0], &u1[0], &_par, &_par2, &out[0]);

    return out;
}

dvec ParBiCop::hFunc2(dvec &u1, dvec &u2)
{
    int n = u1.size();
    dvec out(n);

    Hfunc2(&_family, &n, &u1[0], &u2[0], &_par, &_par2, &out[0]);

    return out;
}

// PDF
dvec ParBiCop::PDF(dvec &u1, dvec &u2)
{
    int n = u1.size();
    dvec out(n);

    PDF_seperate(&_family, &n, &u1[0], &u2[0], &_par, &_par2, &out[0]);

    return out;
};

// fit statistics
double ParBiCop::logLik(dvec &u1, dvec &u2)
{
    int n = u1.size();
    double out;

    LL_mod2(&_family, &n, &u1[0], &u2[0], &_par, &_par2, &out);

    return out;
}

double ParBiCop::AIC(dvec &u1, dvec &u2)
{
    int n = u1.size();
    double out;

    LL_mod2(&_family, &n, &u1[0], &u2[0], &_par, &_par2, &out);
    out = -2 * out + 2 * this->_npars;

    return out;
}

double ParBiCop::BIC(dvec &u1, dvec &u2)
{
    int n = u1.size();
    double out;

    LL_mod2(&_family, &n, &u1[0], &u2[0], &_par, &_par2, &out);
    out = -2 * out + log(n) * this->_npars;

    return out;
}