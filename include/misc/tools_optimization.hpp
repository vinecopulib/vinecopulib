// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <nlopt.hpp>
#include "misc/tools_eigen.hpp"
#include "bicop/bicop_parametric.hpp"

namespace tools_optimization {

    //! A helper struct for nlopt (profile) maximum likelihood estimation
    //!
    typedef struct
    {
        const Eigen::MatrixXd& U; //! The data
        vinecopulib::ParBicop* bicop; //! A pointer to the bivariate copula to optimize
        double par0;  //! The main dependence parameter
        unsigned int objective_calls; //! The number of evaluations of the objective
    } ParBicopOptData;

    //! A class for the controls to nlopt
    //!
    class NLoptControls
    {
    public:

        //! Create controls using the default contructor
        //!
        //! @return The default NLopt controls.
        NLoptControls();

        //! Create controls by passing the arguments
        //!
        //! @param xtol_rel Relative tolerance on the parameter value
        //! @param xtol_abs Absolute tolerance on the parameter value
        //! @param ftol_rel Relative tolerance on the function value
        //! @param ftol_abs Absolue tolerance on the function value
        //! @param maxeval Maximal number of evaluations of the objective
        //! @return Custom NLopt controls.
        NLoptControls(double xtol_rel, double xtol_abs, double ftol_rel, double ftol_abs, int maxeval);


        //! Set controls of an optimizer
        //!
        //! @param opt Pointer to the optimizer for control setting
        void set_controls(nlopt::opt* opt);

        //! Getters and setters.
        //! @{
        double get_xtol_rel();
        double get_xtol_abs();
        double get_ftol_rel();
        double get_ftol_abs();
        double get_maxeval();
        //! @}
    private:
        double xtol_rel_; //! Relative tolerance on the parameter value
        double xtol_abs_; //! Absolute tolerance on the parameter value
        double ftol_rel_; //! Relative tolerance on the function value
        double ftol_abs_; //! Absolute tolerance on the function value
        int maxeval_; //! Maximal number of evaluations of the objective

        //! Sanity checks
        //! @{
        void check_parameters(double xtol_rel, double xtol_abs, double ftol_rel, double ftol_abs, int maxeval);
        //! @}
    };


    class Optimizer {
    public:
        //! Create an optimizer using the default controls
        //!
        //! @param n_parameters Number of parameters to optimize
        //! @return An optimizer using the default controls.
        Optimizer(unsigned int n_parameters);

        //! Create an optimizer using the custom controls
        //!
        //! @param n_parameters Number of parameters to optimize
        //! @param xtol_rel Relative tolerance on the parameter value
        //! @param xtol_abs Absolute tolerance on the parameter value
        //! @param ftol_rel Relative tolerance on the function value
        //! @param ftol_abs Absolue tolerance on the function value
        //! @param maxeval Maximal number of evaluations of the objective
        //! @return An optimizer using the custom controls.
        Optimizer(unsigned int n_parameters, double xtol_rel, double xtol_abs,
                  double ftol_rel, double ftol_abs, int maxeval);

        //! Set the optimizer's bounds
        //!
        //! @param bounds A matrix of parameters bounds
        void set_bounds(
            const Eigen::MatrixXd& lower_bounds,
            const Eigen::MatrixXd& upper_bounds
        );

        //! Set the optimizer's objective and data
        //!
        //! @param f The optimizer's objective function (see nlopt's documentation)
        //! @param data The optimizer's data (see nlopt's documentation)
        void set_objective(nlopt::vfunc f, void* f_data);

        //! Solve the optimization problem
        //!
        //! @param initial_parameters Eigen::Vector of starting values
        //! @return MLE or PMLE
        Eigen::VectorXd optimize(Eigen::VectorXd initial_parameters);

    private:
        unsigned int n_parameters_;
        nlopt::opt opt_;
        NLoptControls controls_;
    };

    // the objective function for maximum likelihood estimation
    double mle_objective(const std::vector<double>& x,
                         std::vector<double>& grad,
                         void* data);

    // the objective function for profile maximum likelihood estimation
    double pmle_objective(const std::vector<double>& x,
                          std::vector<double> &,
                          void* data);

}
