/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VINECOPULIB_BICOP_CLASS_HPP
#define VINECOPULIB_BICOP_CLASS_HPP

#include <Eigen/Dense>
#include <random>
#include <nlopt.hpp>
#include <iostream>
#include <memory>
#include <vector>

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
    //! Create a bivariate copula using the default contructor
    //!
    //! @param family the copula family.
    //! @rotation rotation the rotation type.
    //! @return A pointer to an object that inherits from \c Bicop.
    static std::shared_ptr<Bicop> create(const int& family,
                                         const int& rotation);

    //! Create a bivariate copula with a specified parameters vector
    //!
    //! @param family the copula family.
    //! @param par the copula parameters (must be compatible with family).
    //! @rotation rotation the rotation type.
    //! @return A pointer to an object that inherits from \c Bicop.
    static std::shared_ptr<Bicop> create(const int& family,
                                         const VecXd& parameters,
                                         const int& rotation);

    //! Select a bivariate copula
    //!
    //! @param data the data to fit the bivariate copula.
    //! @param selection_criterion the selection criterion ("aic" or "bic").
    //! @param family_set the set of copula families to consider (if empty, then all families are included).
    //! @param use_rotations whether rotations in the familyset are included.
    //! @param preselect_families whether to exclude families before fitting based on symmetry properties of the data.
    //! @param method indicates the estimation method: either maximum likelihood estimation (method = "mle") or
    //! inversion of Kendall's tau (method = "itau"). When method = "itau" is used with families having more than
    //! one parameter, the main dependence parameter is found by inverting the Kendall's tau and the remainders by a
    //! profile likelihood optimization.
    //! @return A pointer to an object that inherits from \c Bicop.
    static std::shared_ptr<Bicop> select(const MatXd& data,
                                         std::string selection_criterion,
                                         std::vector<int> family_set,
                                         bool use_rotations,
                                         bool preselect_families,
                                         std::string method);

    //! \defgroup df Copula density
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @{
    VecXd pdf(const MatXd& u);
    //! @}

    //! \defgroup hfunctions h-functions
    //!
    //! h-functions are defined as one-dimensional integrals over a bivariate
    //! copula density \f$ c \f$,
    //! \f[ h_1(u_1, u_2) = \int_0^{u_2} c(u_1, s) ds, \f]
    //! \f[ h_2(u_1, u_2) = \int_0^{u_1} c(s, u_2) ds. \f]
    //! \c hinv1 is the inverse w.r.t. second argument (conditioned on first),
    //! \c hinv2 is the inverse w.r.t. first argument (conditioned on second),
    //!
    //! \c hfunc1, \c hfunc2, \c hinv1, and \c hinv2 mainly take care that
    //! rotations are properly handled.  They call \c hfunc1_default,
    //! \c hfunc2_default, \c hinv1_default, and hinv2_default which are
    //! family-specific implementations for `rotation = 0`.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @{
    VecXd hfunc1(const MatXd& u);
    VecXd hfunc2(const MatXd& u);
    VecXd hinv1(const MatXd& u);
    VecXd hinv2(const MatXd& u);
    //! @}

    //! Simulate from a bivariate copula
    //!
    //! @param n number of observations.
    MatXd simulate(const int& n);

    //! \devgroup fitstats Fit methods and statistics
    //!
    //! @param u \f$m \times 2\f$ matrix of observations.
    //! @{
    virtual void fit(const MatXd &data, std::string method) = 0;
    double loglik(const MatXd& u);
    double aic(const MatXd& u);
    double bic(const MatXd& u);
    //! @}

    //! Get number of parameters.
    virtual int calculate_npars() = 0;

    //! Calculate the theoretical Kendall's tau
    double calculate_tau();  // this will be a generic fall back method
    virtual double parameters_to_tau(const VecXd& parameters) = 0;

    //! Getters and setters.
    //! @{
    int get_family() const {return family_;}
    int get_rotation() const {return rotation_;}
    std::string get_association_direction() const {return association_direction_;}
    VecXd get_parameters() const {return parameters_;}
    MatXd get_parameters_bounds() const {return parameter_bounds_;}

    void set_rotation(const int& rotation);
    void set_parameters(const VecXd& parameters) {parameters_ = parameters;}
    //! @}


protected:
    virtual VecXd pdf_default(const MatXd& u) = 0;
    virtual VecXd hfunc1_default(const MatXd& u) = 0;
    virtual VecXd hfunc2_default(const MatXd& u) = 0;
    virtual VecXd hinv1_default(const MatXd& u) = 0;
    virtual VecXd hinv2_default(const MatXd& u) = 0;

    //! Numerical inversion of h-functions
    //!
    //! These are generic functions to invert the hfunctions numerically.
    //! They can be used in derived classes to define \c hinv1 and \c hinv2.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    VecXd hinv1_num(const MatXd& u);
    VecXd hinv2_num(const MatXd& u);

    //! Data manipulations for rotated families
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    MatXd rotate_u(const MatXd& u);
    MatXd swap_cols(const MatXd& u);

    int family_;
    int rotation_;
    std::string association_direction_;
    VecXd parameters_;
    MatXd parameter_bounds_;   // first row lower, second row upper
};

typedef std::shared_ptr<Bicop> BicopPtr;
double correlation(const MatXd& z);
template<typename T> bool is_member(T element, std::vector<T> set);
std::vector<double> get_c1c2(const MatXd& data, double tau);
bool preselect_family(double c1, double c2, double tau, int family, int rotation, bool is_rotationless);

#endif
