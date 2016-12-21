//
// Created by Vatter Thibault on 21/12/16.
//

#ifndef VINECOPLIB_BICOP_STUDENT_HPP
#define VINECOPLIB_BICOP_STUDENT_HPP

#include "bicop_elliptical.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

class StudentBicop : public EllipticalBicop {

public:
    // constructor
    StudentBicop();
    StudentBicop(double rho, double nu);

    // bounds on the copula parameter
    MatXd get_bounds();

    // hfunction
    VecXd hfunc1(const MatXd &u);

    // inverse hfunction
    VecXd hinv1(const MatXd &u);

    // PDF
    VecXd pdf(const MatXd &u);
};


#endif //VINECOPLIB_BICOP_STUDENT_HPP
