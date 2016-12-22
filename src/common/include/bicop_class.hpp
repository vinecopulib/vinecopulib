/*
Copyright 2016 Thibault Vatter

This file is part of vinecopulib.

vinecoplib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecoplib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VINECOPLIB_BICOP_CLASS_HPP
#define VINECOPLIB_BICOP_CLASS_HPP

#include <Eigen/Dense>
#include <random>
#include <nlopt.hpp>
#include <iostream>
typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

extern "C" {
    #include "ktau.h"
}

//! A class for bivariate copulas
//!
//! The Bicop class is abstract, you cannot create an instance of this class,
//! but only of the derived classes.
class Bicop
{

public:
    //! Create a bivariate copula
    //!
    //! @param family the copula family.
    //! @param par the copula parameters (must be compatible with family).
    //! @return A pointer to an object that inherits from \c Bicop.
    static Bicop* create(const int& family, const VecXd& parameters);

    //! Copula density
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    virtual VecXd pdf(const MatXd& u) = 0;

    //! \defgroup hfunctions h-functions
    //!
    //! h-functions are defined as one-dimensional integrals over a bivariate
    //! copula density \f$ c \f$,
    //! \f[ h_1(u_1, u_2) = \int_0^{u_2} c(u_1, s) ds, \f]
    //! \f[ h_2(u_1, u_2) = \int_0^{u_1} c(s, u_2) ds. \f]
    //! \c hinv1 is the inverse w.r.t. second argument (conditioned on first),
    //! \c hinv2 is the inverse w.r.t. first argument (conditioned on second),
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @{
    virtual VecXd hfunc1(const MatXd& u) = 0;
    virtual VecXd hfunc2(const MatXd& u) = 0;
    virtual VecXd hinv1(const MatXd& u) = 0;
    virtual VecXd hinv2(const MatXd& u) = 0;
    //! @}

    //! Numerical inversion of h-functions
    //!
    //! These are generic functions to invert the hfunctions numerically.
    //! They can be used in derived classes to define \c hinv1 and \c hinv2.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    VecXd hinv1_num(const MatXd& u);
    VecXd hinv2_num(const MatXd& u);

    //! Simulate from a bivariate copula
    //!
    //! @param n number of observations.
    MatXd simulate(const int& n);

    //! \devgroup fitstats Fit statistics
    //!
    //! @param u \f$m \times 2\f$ matrix of observations.
    //! @{
    double loglik(MatXd& u);
    double aic(MatXd& u);
    double bic(MatXd& u);
    //! @}

    //! Get number of parameters.
    virtual int calculate_npars() = 0;

    //! Calculate the theoretical Kendall's tau
    double calculate_tau();  // this will be a generic fall back method
    virtual double par_to_tau(const VecXd& parameters);

    //! Getters and setters.
    //! @{
    int get_family() const {return family_;}
    int get_rotation() const {return rotation_;}
    VecXd get_parameters() const {return parameters_;}
    VecXd get_par_bounds() const {return parameter_bounds_;}

    void set_rotation(const int& rotation) {rotation_ = rotation;}
    void set_parameters(const VecXd& parameters) {parameters_ = parameters;}
    //! @}

protected:
    int family_;
    int rotation_;
    VecXd parameters_;
    MatXd parameter_bounds_;
};

#endif
