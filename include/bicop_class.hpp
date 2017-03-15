// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <memory>
#include <vector>

#include "tools_eigen.hpp"

namespace vinecopulib
{
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
                                             const VectorXd& parameters,
                                             const int& rotation);

        static std::shared_ptr<Bicop> select(const MatrixXd& data,
                                             std::vector<int> family_set = {0, 1, 2, 3, 4, 5, 6, 1001},
                                             std::string method = "mle",
                                             std::string selection_criterion = "bic",
                                             bool preselect_families = true);

        VectorXd pdf(const MatrixXd& u);
        VectorXd hfunc1(const MatrixXd& u);
        VectorXd hfunc2(const MatrixXd& u);
        VectorXd hinv1(const MatrixXd& u);
        VectorXd hinv2(const MatrixXd& u);
        MatrixXd simulate(const int& n);

        virtual void fit(const MatrixXd &data, std::string method) = 0;
        double loglik(const MatrixXd& u);
        double aic(const MatrixXd& u);
        double bic(const MatrixXd& u);

        virtual double calculate_npars() = 0;
        double calculate_tau();  // this will be a generic fall back method
        virtual double parameters_to_tau(const VectorXd& parameters) = 0;
        VectorXd tau_to_parameters(const double& tau);


        int get_family() const;
        int get_rotation() const;
        std::string get_association_direction() const;
        VectorXd get_parameters() const;
        MatrixXd get_parameters_bounds() const;
        void set_rotation(const int& rotation);
        void set_parameters(const VectorXd& parameters);


        //! Adjust the copula to a flipping of arguments (u,v) -> (v,u)
        virtual void flip() = 0;

    protected:
        virtual VectorXd pdf_default(const MatrixXd& u) = 0;
        virtual VectorXd hfunc1_default(const MatrixXd& u) = 0;
        virtual VectorXd hfunc2_default(const MatrixXd& u) = 0;
        virtual VectorXd hinv1_default(const MatrixXd& u) = 0;
        virtual VectorXd hinv2_default(const MatrixXd& u) = 0;

        VectorXd hinv1_num(const MatrixXd& u);
        VectorXd hinv2_num(const MatrixXd& u);

        MatrixXd cut_and_rotate(const MatrixXd& u);
        MatrixXd swap_cols(const MatrixXd& u);

        int family_;
        std::string family_name_;
        int rotation_;
        std::string association_direction_;
        VectorXd parameters_;
        MatrixXd parameters_bounds_;   // first column lower, second column upper

    private:
        void check_parameters(const VectorXd& parameters);
        void check_rotation(const int& rotation);
    };

    typedef std::shared_ptr<Bicop> BicopPtr;
}
