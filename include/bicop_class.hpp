/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include "distribution_tools.hpp"
#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <memory>
#include <vector>
#include <exception>
#include <functional>

#define _USE_MATH_DEFINES
#include <cmath>

typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

extern "C" {
#include "c_tools.h"
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
    //! @param rotation the rotation type.
    //! @return A pointer to an object that inherits from \c Bicop.
    static std::shared_ptr<Bicop> create(const int& family,
                                         const int& rotation);

    //! Create a bivariate copula with a specified parameters vector
    //!
    //! @param family the copula family.
    //! @param par the copula parameters (must be compatible with family).
    //! @param rotation the rotation type.
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
    virtual double calculate_npars() = 0;

    //! Calculate the theoretical Kendall's tau
    double calculate_tau();  // this will be a generic fall back method
    virtual double parameters_to_tau(const VecXd& parameters) = 0;

    //! Extract the copula parameter from Kendall's tau whenever possible
    VecXd tau_to_parameters(const double& tau);

    //! Getters and setters.
    //! @{
    int get_family() const {return family_;}
    int get_rotation() const {return rotation_;}
    std::string get_association_direction() const {return association_direction_;}
    VecXd get_parameters() const {return parameters_;}
    MatXd get_parameters_bounds() const {return parameters_bounds_;}

    void set_rotation(const int& rotation);
    void set_parameters(const VecXd& parameters);
    //! @}

    //! Adjust the copula to a flipping of arguments (u,v) -> (v,u)
    virtual void flip() = 0;

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
    MatXd cut_and_rotate(const MatXd& u);
    MatXd swap_cols(const MatXd& u);

    int family_;
    std::string family_name_;
    int rotation_;
    std::string association_direction_;
    VecXd parameters_;
    MatXd parameters_bounds_;   // first column lower, second column upper

private:
    //! Sanity checks
    //! @{
    void check_parameters(const VecXd& parameters);
    void check_rotation(const int& rotation);
    //! @}
};

typedef std::shared_ptr<Bicop> BicopPtr;
double correlation(const MatXd& z);
std::vector<double> get_c1c2(const MatXd& data, double tau);
bool preselect_family(double c1, double c2, double tau, int family, int rotation, bool is_rotationless);
VecXd invert_f(const VecXd &x, std::function<VecXd(const VecXd&)> f, const double lb = 1e-20, const double ub = 1-1e-20,
               int n_iter = 35);

template<typename T> bool is_member(T element, std::vector<T> set) {
    return std::find(set.begin(), set.end(), element) != set.end(); }


template<typename T> std::vector<T> intersect(std::vector<T> x, std::vector<T> y) {
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    std::vector<T> common;
    std::set_intersection(
            x.begin(), x.end(),
            y.begin(), y.end(),
            std::back_inserter(common)
    );
    return common;
}
