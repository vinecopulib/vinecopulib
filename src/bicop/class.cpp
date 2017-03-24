// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/class.hpp"
#include "bicop/tools_bicopselect.hpp"
#include "misc/tools_stats.hpp"
#include "misc/tools_stl.hpp"

#include <iostream>

namespace vinecopulib
{
    
    /// @name Constructors
    /// 
    /// There are several ways to construct a `Bicop` object. The copula model
    /// is fully characterized by the family, rotation, and parameters. You can
    /// construct a model directly from these characteristics or use the 
    /// data-driven constructor (see select()).
    /// 
    /// @param family the copula family.
    /// @param rotation the rotation of the copula; one of 0, 90, 180, or 270 
    ///     (for Independence, Gaussian, Student, Frank, and nonparametric 
    ///     families, only 0 is allowed).
    /// @param parameters the copula parameters.
    /// @param data see select().
    /// @param family_set see select().
    /// @param method see select().
    /// @param selection_criterion see select().
    /// @param preselect_families see select().
    /// 
    /// @{
    Bicop::Bicop()
    {
        bicop_ = AbstractBicop::create();
    }
    
    Bicop::Bicop(BicopFamily family, int rotation, 
        const Eigen::MatrixXd& parameters)
    {
        bicop_ = AbstractBicop::create(family, parameters);
        // family must be set before checking the rotation
        set_rotation(rotation);
    }
    
    Bicop::Bicop(
            Eigen::Matrix<double, Eigen::Dynamic, 2> data,
            std::vector<BicopFamily> family_set,
            std::string method,
            std::string selection_criterion,
            bool preselect_families)
    {
        select(data, family_set, method, selection_criterion, preselect_families);
    }


    //! evaluates the copula density.
    //!
    //! @param u \f$n \times 2\f$ matrix of evaluation points.
    //! @return The copula density evaluated at \c u.
    Eigen::VectorXd Bicop::pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        Eigen::VectorXd f = bicop_->pdf(cut_and_rotate(u));
        f = f.unaryExpr([](const double x){ return std::min(x,1e16);});
        return f;
    }

    //! @name h-functions
    //!
    //! calculates h-functions and their inverse. h-functions are defined as 
    //! one-dimensional integrals over a bivariate copula density \f$ c \f$,
    /// \f{align*}{
    ///     h_1(u_1, u_2) = \int_0^{u_2} c(u_1, s) ds, \qquad
    ///     h_2(u_1, u_2) = \int_0^{u_1} c(s, u_2) ds.
    ///  \f}
    //! `hinv1` is the inverse w.r.t. the second argument (conditioned on first),
    //! `hinv2` is the inverse w.r.t. the first argument (conditioned on second).
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The (inverse) h-function evaluated at `u`.
    //! @{
    Eigen::VectorXd Bicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (rotation_) {
            case 0:
                return bicop_->hfunc1(cut_and_rotate(u));

            case 90:
                return bicop_->hfunc2(cut_and_rotate(u));

            case 180:
                return 1.0 - bicop_->hfunc1(cut_and_rotate(u)).array();

            case 270:
                return 1.0 - bicop_->hfunc2(cut_and_rotate(u)).array();

            default:
                throw std::runtime_error(std::string(
                        "rotation can only take values in {0, 90, 180, 270}"
                ));
        }
    }

    Eigen::VectorXd Bicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (rotation_) {
            case 0:
                return bicop_->hfunc2(cut_and_rotate(u));

            case 90:
                return 1.0 - bicop_->hfunc1(cut_and_rotate(u)).array();

            case 180:
                return 1.0 - bicop_->hfunc2(cut_and_rotate(u)).array();

            case 270:
                return bicop_->hfunc1(cut_and_rotate(u));

            default:
                throw std::runtime_error(std::string(
                        "rotation can only take values in {0, 90, 180, 270}"
                ));
        }
    }

    Eigen::VectorXd Bicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (rotation_) {
            case 0:
                return bicop_->hinv1(cut_and_rotate(u));

            case 90:
                return bicop_->hinv2(cut_and_rotate(u));

            case 180:
                return 1.0 - bicop_->hinv1(cut_and_rotate(u)).array();

            case 270:
                return 1.0 - bicop_->hinv2(cut_and_rotate(u)).array();

            default:
                throw std::runtime_error(std::string(
                        "rotation only takes value in {0, 90, 180, 270}"
                ));
        }
    }

    Eigen::VectorXd Bicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (rotation_) {
            case 0:
                return bicop_->hinv2(cut_and_rotate(u));

            case 90:
                return 1.0 - bicop_->hinv1(cut_and_rotate(u)).array();

            case 180:
                return 1.0 - bicop_->hinv2(cut_and_rotate(u)).array();

            case 270:
                return bicop_->hinv1(cut_and_rotate(u));

            default:
                throw std::runtime_error(std::string(
                        "rotation only takes value in {0, 90, 180, 270}"
                ));
        }
    }
    //! @}


    //! simulates from a bivariate copula.
    //!
    //! @param n number of observations.
    //! @return An \f$ n \times 2 \f$ matrix of samples from the copula model.
    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::simulate(const int& n)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> U =
                tools_stats::simulate_uniform(n, 2);
        // use inverse Rosenblatt transform to generate a sample from the copula
        U.col(1) = hinv1(U);
        return U;
    }

    //! @anchor fit_statistics
    //! @name Fit statistics
    //! 
    //! Calculates fit statistics for a bivariate copula. The log-likelihood,
    //! Aikaike information criterion (AIC), and Bayesian information criterion
    //! (BIC) are defined as
    /// \f{align*}{
    ///     \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, U_{2, i}), \quad
    ///     \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \quad
    ///     \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p,     
    /// \f}
    /// where \f$ c \f$ is the copula density pdf() and  \f$ p \f$ is the 
    /// (effective) number of parameters of the model. 
    /// The loglikehood is used for model fitting, the others for selecting 
    /// between fitted models.  The AIC is consistent for nonparametric models,
    /// the BIC is more conservative and consistent for parametric models.
    /// The effective number of parameters is the actual number of parameters 
    //! for parameteric families. For nonparametric families, there is another
    //! definition that is conceptually similar in the sense that it can be used
    //! in the calculation of fit statistics.
    //! 
    //! 
    //! @param u \f$n \times 2\f$ matrix of observations.
    //! @return The fit statistic for the model and observations `u`.
    //! @{
    double Bicop::loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return pdf(u).array().log().sum();
    }

    double Bicop::aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * loglik(u) + 2 * calculate_npars();
    }

    double Bicop::bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * loglik(u) + calculate_npars() * log(u.rows());
    }
    
    double Bicop::calculate_npars()
    {
        return bicop_->calculate_npars();
    }
    //! @}
    
    
    //! @name Conversion between Parameters and Kendall's tau 
    //! 
    //! Currently `parameters_to_tau()` works for all families but 
    //! `BicopFamily::tll0`. `tau_to_parameters()` only works for one-parameter
    //! families.
    //! 
    //! @param tau 
    //! @param parameters 
    //! @return The tau/parameters corresponding to parameters/tau for the 
    //!     current family.
    //! @{
    Eigen::MatrixXd Bicop::tau_to_parameters(const double& tau)
    {
        return bicop_->tau_to_parameters(tau);
    }

    double Bicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        double tau = bicop_->parameters_to_tau(parameters);
        if (tools_stl::is_member(rotation_, {90, 270})) {
            tau *= -1;
        }
        return tau;
    }
    //! @}
    
    //! @name Getters and setters
    //! 
    //! @param rotation 
    //! @param parameters 
    //! 
    //! @{
    BicopFamily Bicop::get_family() const 
    {
        return bicop_->get_family();
    }
    
    std::string Bicop::get_family_name() const 
    {
        return bicop_->get_family_name();
    };
    
    int Bicop::get_rotation() const
    {
        return rotation_;
    }
    
    Eigen::MatrixXd Bicop::get_parameters() const 
    {
        return bicop_->get_parameters();
    }

    void Bicop::set_rotation(int rotation) {
        check_rotation(rotation);
        rotation_ = rotation;
    }

    void Bicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        bicop_->set_parameters(parameters);
    }
    //! @}
    
    
    //! @name Utilities
    //! 
    //! `flip()` adjust's the copula model to a change in the variable order.
    //! `str()` summarizes the model into a string (can be used for printing).
    void Bicop::flip()
    {
        BicopFamily family = bicop_->get_family();
        if (tools_stl::is_member(family, bicop_families::flip_by_rotation)) {
            if (rotation_ == 90) {
                set_rotation(270);
            } else if (rotation_ == 270) {
                set_rotation(90);
            }
        } else {
            bicop_->flip();    
        }
    }
    
    std::string Bicop::str()
    {
        std::stringstream bicop_str;
        bicop_str << "family = "    << get_family_name() <<
                  ", rotation = "   << get_rotation() <<
                  ", parameters = " << get_parameters();

        return bicop_str.str().c_str();
    }
    //! @}

    BicopPtr Bicop::get_bicop()
    {
        return bicop_;
    };
    
    
    //! Estimation and selection of bivariate copula models
    //! 
    //! For parametric models, two different methods are available. `"mle"` fits
    //! the parameters by maximum-likelihood. `"itau"` uses inversion of 
    //! Kendall's \f$ \tau \f$, but is only available for one-parameter families
    //! and the Student t copula. For the latter, there is a one-to-one 
    //! transformation for the first parameter, the second is found by profile
    //! likelihood optimization (with accuracy of at least 0.5). nonparametric
    //! families have specialized methods, no specification is required.
    //! 
    //! @param data an \f$ n \times 2 \f$ matrix of observations contained in
    //!     \f$(0, 1)^2 \f$.
    //! @param family_set the set of copula families to consider (if empty, then
    //!     all families are included).
    //! @param method the fit method; possible choices: `"mle"`, `"itau"`.
    //! @param selection_criterion the selection criterion (`"aic"` or `"bic"`).
    //! @param preselect_families whether to exclude families before fitting
    //!     based on symmetry properties of the data.
    //! @{
    void Bicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
            std::string method)
    {
        bicop_->fit(cut_and_rotate(data), method);
    }

    void Bicop::select(
            Eigen::Matrix<double, Eigen::Dynamic, 2> data,
            std::vector<BicopFamily> family_set,
            std::string method,
            std::string selection_criterion,
            bool preselect_families
    )
    {
        using namespace tools_stl;
        // If the familyset is empty, use all families.
        // If the familyset is not empty, check that all included families are implemented.
        if (family_set.empty()) {
            if (method == "itau") {
                family_set = bicop_families::itau;
            } else {
                family_set = bicop_families::all;
            }
        } else {
            if (intersect(family_set, bicop_families::all).empty()) {
                throw std::runtime_error(
                        std::string("One of the families is not implemented")
                );
            }
            if (method == "itau") {
                family_set = intersect(family_set, bicop_families::itau);
                if (family_set.empty()) {
                    throw std::runtime_error(
                            std::string("No family with method itau provided")
                    );
                }
            }
        }

        // When using rotations, add only the ones that yield the appropriate
        // association direction.
        auto tau = tools_stats::pairwise_ktau(data);
        std::vector<int> which_rotations;
        if (tau > 0) {
            which_rotations = {0, 180};
        } else {
            which_rotations = {90, 270};
        }

        std::vector<double> c(2);
        if (preselect_families) {
            c = get_c1c2(data, tau);
        }

        // Create the combinations of families and rotations to estimate
        std::vector<BicopFamily> families;
        std::vector<int> rotations;
        for (auto family : family_set) {
            bool is_rotationless = is_member(family,
                                             bicop_families::rotationless);
            bool preselect = true;
            if (is_rotationless) {
                if (preselect_families) {
                    preselect = preselect_family(c, tau, family,
                                                 0, is_rotationless);
                }
                if (preselect) {
                    families.push_back(family);
                    rotations.push_back(0);
                }
            } else {
                for (auto rotation : which_rotations) {
                    if (preselect_families) {
                        preselect = preselect_family(c, tau, family,
                                                     rotation, is_rotationless);
                    }
                    if (preselect) {
                        families.push_back(family);
                        rotations.push_back(rotation);
                    }
                }
            }
        }

        // Estimate all models and select the best one using the selection_criterion
        BicopPtr fitted_bicop;
        int fitted_rotation;
        double fitted_criterion = 1e6;
        for (unsigned int j = 0; j < families.size(); j++) {
            // Estimate the model
            bicop_ = AbstractBicop::create(families[j]);
            rotation_ = rotations[j];
            bicop_->fit(cut_and_rotate(data), method);

            // Compute the selection criterion
            double new_criterion;
            if (selection_criterion == "aic") {
                new_criterion = aic(data);
            } else if (selection_criterion == "bic") {
                new_criterion = bic(data);
            } else {
                throw std::runtime_error("Selection criterion not implemented");
            }

            // If the new model is better than the current one, then replace the current model by the new one
            if (new_criterion < fitted_criterion) {
                fitted_criterion = new_criterion;
                fitted_bicop = bicop_;
                fitted_rotation = rotation_;
            }
        }
        
        bicop_ = fitted_bicop;
        rotation_ = fitted_rotation;
    }
    //! @}
    
    //! Data manipulations for rotated families
    //!
    //! @param u \f$m \times 2\f$ matrix of data.
    //! @return The manipulated data.
    //! @{
    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::cut_and_rotate(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> u_new(u.rows(), 2);

        // counter-clockwise rotations
        switch (rotation_) {
            case 0:
                u_new = u;
                break;

            case 90:
                u_new.col(0) = u.col(1);
                u_new.col(1) = 1.0 - u.col(0).array();
                break;

            case 180:
                u_new.col(0) = 1.0 - u.col(0).array();
                u_new.col(1) = 1.0 - u.col(1).array();
                break;

            case 270:
                u_new.col(0) = 1.0 - u.col(1).array();
                u_new.col(1) = u.col(0);
                break;
        }

        // truncate to interval [eps, 1 - eps]
        Eigen::Matrix<double, Eigen::Dynamic, 2> eps = 
            Eigen::Matrix<double, Eigen::Dynamic, 2>::Constant(u.rows(), 2, 1e-10);
        u_new = u_new.array().min(1.0 - eps.array());
        u_new = u_new.array().max(eps.array());

        return u_new;
    }
    //! @}
    
    void Bicop::check_rotation(int rotation)
    {
        using namespace tools_stl;
        std::vector<int> allowed_rotations = {0, 90, 180, 270};
        if (!is_member(rotation, allowed_rotations)) {
            throw std::runtime_error("rotation must be one of {0, 90, 180, 270}");
        }
        if (is_member(bicop_->get_family(), bicop_families::rotationless)) {
            if (rotation != 0) {
                throw std::runtime_error("rotation must be 0 for the " + 
                    bicop_->get_family_name() + " copula");
            }
        }
    }
}
