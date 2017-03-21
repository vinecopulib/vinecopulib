// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <memory>
#include <vector>
#include "tools_eigen.hpp"
#include "bicop_family.hpp"

namespace vinecopulib
{
    //! A class for bivariate copulas
    //!
    //! The Bicop class is abstract, you cannot create an instance of this class,
    //! but only of the derived classes.
    class Bicop
    {
    public:
        static std::shared_ptr<Bicop> create(
            BicopFamily family = BicopFamily::indep,
            int rotation = 0
        );
        
        static std::shared_ptr<Bicop> create(
            BicopFamily family,
            int rotation,
            Eigen::VectorXd parameters
        );

        static std::shared_ptr<Bicop> select(
            Eigen::MatrixXd& data,
            std::vector<BicopFamily> family_set = bicop_families::all,
            std::string method = "mle",
            std::string selection_criterion = "bic",
            bool preselect_families = true
        );

        Eigen::VectorXd pdf(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u);
        Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv1(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2(const Eigen::MatrixXd& u);
        Eigen::MatrixXd simulate(const int& n);

        virtual void fit(const Eigen::MatrixXd &data, std::string method) = 0;
        double loglik(const Eigen::MatrixXd& u);
        double aic(const Eigen::MatrixXd& u);
        double bic(const Eigen::MatrixXd& u);

        virtual double calculate_npars() = 0;
        double calculate_tau();  // this will be a generic fall back method
        virtual double parameters_to_tau(const Eigen::VectorXd& parameters) = 0;
        Eigen::VectorXd tau_to_parameters(const double& tau);


        BicopFamily get_family() const;
        std::string get_family_name() const;
        int get_rotation() const;
        std::string get_association_direction() const;
        Eigen::VectorXd get_parameters() const;
        Eigen::MatrixXd get_parameters_bounds() const;
        void set_rotation(const int& rotation);
        void set_parameters(const Eigen::VectorXd& parameters);


        //! Adjust the copula to a flipping of arguments (u,v) -> (v,u)
        virtual void flip() = 0;

    protected:
        virtual Eigen::VectorXd pdf_default(const Eigen::MatrixXd& u) = 0;
        virtual Eigen::VectorXd hfunc1_default(const Eigen::MatrixXd& u) = 0;
        virtual Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u) = 0;
        virtual Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u) = 0;
        virtual Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u) = 0;
        virtual Eigen::VectorXd tau_to_parameters_default(const double& tau) = 0;

        Eigen::VectorXd hinv1_num(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2_num(const Eigen::MatrixXd& u);

        Eigen::MatrixXd cut_and_rotate(const Eigen::MatrixXd& u);
        Eigen::MatrixXd swap_cols(const Eigen::MatrixXd& u);

        BicopFamily family_;
        int rotation_;
        Eigen::VectorXd parameters_;
        Eigen::MatrixXd parameters_bounds_;   // first column lower, second column upper

    private:
        void check_parameters(const Eigen::VectorXd& parameters);
        void check_rotation(const int& rotation);
    };

    typedef std::shared_ptr<Bicop> BicopPtr;
    Eigen::VectorXd no_tau_to_parameters(const double&);
}
