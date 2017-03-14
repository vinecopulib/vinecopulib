// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <memory>
#include <vector>
#include "tools_eigen.hpp"


//! A class for bivariate copulas
//!
//! The Bicop class is abstract, you cannot create an instance of this class,
//! but only of the derived classes.
class Bicop
{

public:
    static std::shared_ptr<Bicop> create(const int& family,
                                         const int& rotation);
    static std::shared_ptr<Bicop> create(const int& family,
                                         const VecXd& parameters,
                                         const int& rotation);

    static std::shared_ptr<Bicop> select(const MatXd& data,
                                         std::vector<int> family_set = {0, 1, 2, 3, 4, 5, 6, 1001},
                                         std::string method = "mle",
                                         std::string selection_criterion = "bic",
                                         bool preselect_families = true);

    VecXd pdf(const MatXd& u);
    VecXd hfunc1(const MatXd& u);
    VecXd hfunc2(const MatXd& u);
    VecXd hinv1(const MatXd& u);
    VecXd hinv2(const MatXd& u);
    MatXd simulate(const int& n);

    virtual void fit(const MatXd &data, std::string method) = 0;
    double loglik(const MatXd& u);
    double aic(const MatXd& u);
    double bic(const MatXd& u);

    virtual double calculate_npars() = 0;
    double calculate_tau();  // this will be a generic fall back method
    virtual double parameters_to_tau(const VecXd& parameters) = 0;
    VecXd tau_to_parameters(const double& tau);


    int get_family() const;
    int get_rotation() const;
    std::string get_association_direction() const;
    VecXd get_parameters() const;
    MatXd get_parameters_bounds() const;
    void set_rotation(const int& rotation);
    void set_parameters(const VecXd& parameters);


    //! Adjust the copula to a flipping of arguments (u,v) -> (v,u)
    virtual void flip() = 0;

protected:
    virtual VecXd pdf_default(const MatXd& u) = 0;
    virtual VecXd hfunc1_default(const MatXd& u) = 0;
    virtual VecXd hfunc2_default(const MatXd& u) = 0;
    virtual VecXd hinv1_default(const MatXd& u) = 0;
    virtual VecXd hinv2_default(const MatXd& u) = 0;

    VecXd hinv1_num(const MatXd& u);
    VecXd hinv2_num(const MatXd& u);

    MatXd cut_and_rotate(const MatXd& u);
    MatXd swap_cols(const MatXd& u);

    int family_;
    std::string family_name_;
    int rotation_;
    std::string association_direction_;
    VecXd parameters_;
    MatXd parameters_bounds_;   // first column lower, second column upper

private:
    void check_parameters(const VecXd& parameters);
    void check_rotation(const int& rotation);
};

typedef std::shared_ptr<Bicop> BicopPtr;
std::vector<double> get_c1c2(const MatXd& data, double tau);
bool preselect_family(double c1, double c2, double tau, int family, int rotation, bool is_rotationless);
