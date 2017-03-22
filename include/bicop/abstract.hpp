// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <memory>
#include <vector>

#include "misc/tools_eigen.hpp"
#include "bicop/family.hpp"

namespace vinecopulib
{
    //! A class for bivariate copulas
    //!
    //! The Bicop class is abstract, you cannot create an instance of this class,
    //! but only of the derived classes.
    class AbstractBicop
    {
    public:
        static std::shared_ptr<AbstractBicop> create(
            BicopFamily family = BicopFamily::indep,
            int rotation = 0
        );
        
        static std::shared_ptr<AbstractBicop> create(
            BicopFamily family,
            int rotation,
            Eigen::VectorXd parameters
        );

        static std::shared_ptr<AbstractBicop> select(
            Eigen::Matrix<double, Eigen::Dynamic, 2> data,
            std::vector<BicopFamily> family_set = bicop_families::all,
            std::string method = "mle",
            std::string selection_criterion = "bic",
            bool preselect_families = true
        );

        Eigen::VectorXd pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::Matrix<double, Eigen::Dynamic, 2> simulate(const int& n);

        virtual void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data, std::string method) = 0;
        double loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        double aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        double bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);

        virtual double calculate_npars() = 0;
        virtual double parameters_to_tau(const Eigen::VectorXd& parameters) = 0;
        Eigen::MatrixXd tau_to_parameters(const double& tau);


        BicopFamily get_family() const;
        std::string get_family_name() const;
        int get_rotation() const;
        std::string get_association_direction() const;
        Eigen::MatrixXd get_parameters() const;
        Eigen::MatrixXd get_parameters_lower_bounds() const;
        Eigen::MatrixXd get_parameters_upper_bounds() const;
        void set_rotation(const int& rotation);
        void set_parameters(const Eigen::MatrixXd& parameters);


        //! Adjust the copula to a flipping of arguments (u,v) -> (v,u)
        virtual void flip() = 0;

    protected:
        virtual Eigen::VectorXd pdf_default(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hfunc1_default(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hfunc2_default(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hinv1_default(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hinv2_default(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::MatrixXd tau_to_parameters_default(const double& tau) = 0;

        Eigen::VectorXd hinv1_num(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hinv2_num(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);

        Eigen::Matrix<double, Eigen::Dynamic, 2> cut_and_rotate(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::Matrix<double, Eigen::Dynamic, 2> swap_cols(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);

        BicopFamily family_;
        int rotation_;
        Eigen::MatrixXd parameters_;
        Eigen::MatrixXd parameters_lower_bounds_;
        Eigen::MatrixXd parameters_upper_bounds_;

    private:
        void check_parameters(const Eigen::MatrixXd& parameters);
        void check_parameters_size(const Eigen::MatrixXd& parameters);
        void check_parameters_upper(const Eigen::MatrixXd& parameters);
        void check_parameters_lower(const Eigen::MatrixXd& parameters);
        void check_rotation(const int& rotation);
    };

    typedef std::shared_ptr<AbstractBicop> BicopPtr;
    Eigen::VectorXd no_tau_to_parameters(const double&);
}
