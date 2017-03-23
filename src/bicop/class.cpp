// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/class.hpp"
#include "bicop/tools_bicopselect.hpp"
#include "misc/tools_stats.hpp"
#include "misc/tools_stl.hpp"

namespace vinecopulib
{
    Bicop::Bicop()
    {
        bicop_ = AbstractBicop::create();
    }

    Bicop::Bicop(BicopFamily family, int rotation)
    {
        bicop_ = AbstractBicop::create(family, rotation);
    }
    
    Bicop::Bicop(BicopFamily family, int rotation, Eigen::VectorXd parameters)
    {
        bicop_ = AbstractBicop::create(family, rotation, parameters);
    }

    //! Select a bivariate copula
    //!
    //! @param data the data to fit the bivariate copula.
    //! @param selection_criterion the selection criterion ("aic" or "bic").
    //! @param family_set the set of copula families to consider (if empty, then
    //!     all families are included).
    //! @param use_rotations whether rotations in the familyset are included.
    //! @param preselect_families whether to exclude families before fitting
    //!     based on symmetry properties of the data.
    //! @param method indicates the estimation method: either maximum likelihood
    //!     estimation (method = "mle") or inversion of Kendall's tau (method =
    //!     "itau"). When method = "itau" is used with families having more than
    //!     one parameter, the main dependence parameter is found by inverting
    //!     the Kendall's tau and the remainders by a profile likelihood
    //!     optimization.
    //! @return A pointer to an object that inherits from \c Bicop.
    Bicop::Bicop(
            Eigen::Matrix<double, Eigen::Dynamic, 2> data,
            std::vector<BicopFamily> family_set,
            std::string method,
            std::string selection_criterion,
            bool preselect_families)
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
        double fitted_criterion = 1e6;
        for (unsigned int j = 0; j < families.size(); j++) {
            // Estimate the model
            bicop_ = AbstractBicop::create(families[j], rotations[j]);
            bicop_->fit(data, method);

            // Compute the selection criterion
            double new_criterion;
            if (selection_criterion == "aic") {
                new_criterion = aic(data);
            } else if (selection_criterion == "bic") {
                new_criterion = bic(data);
            } else {
                throw std::runtime_error(std::string("Selection criterion not implemented"));
            }

            // If the new model is better than the current one, then replace the current model by the new one
            if (new_criterion < fitted_criterion) {
                fitted_criterion = new_criterion;
                fitted_bicop = bicop_;
            }
        }
        bicop_ = fitted_bicop;
    }


    //! Copula density
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The copula density evaluated at \c u.
    Eigen::VectorXd Bicop::pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        Eigen::VectorXd f = bicop_->pdf(bicop_->cut_and_rotate(u));
        f = f.unaryExpr([](const double x){ return std::min(x,1e16);});
        return f;
    }

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
    //! rotations are properly handled.  They call \c hfunc1,
    //! \c hfunc2, \c hinv1_default, and hinv2_default which are
    //! family-specific implementations for `rotation = 0`.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The (inverse) h-function evaluated at \c u.
    //! @{
    Eigen::VectorXd Bicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (bicop_->rotation_) {
            case 0:
                return bicop_->hfunc1(bicop_->cut_and_rotate(u));

            case 90:
                return bicop_->hfunc2(bicop_->cut_and_rotate(u));

            case 180:
                return 1.0 - bicop_->hfunc1(bicop_->cut_and_rotate(u)).array();

            case 270:
                return 1.0 - bicop_->hfunc2(bicop_->cut_and_rotate(u)).array();

            default:
                throw std::runtime_error(std::string(
                        "rotation can only take values in {0, 90, 180, 270}"
                ));
        }
    }

    Eigen::VectorXd Bicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (bicop_->rotation_) {
            case 0:
                return bicop_->hfunc2(bicop_->cut_and_rotate(u));

            case 90:
                return 1.0 - bicop_->hfunc1(bicop_->cut_and_rotate(u)).array();

            case 180:
                return 1.0 - bicop_->hfunc2(bicop_->cut_and_rotate(u)).array();

            case 270:
                return bicop_->hfunc1(bicop_->cut_and_rotate(u));

            default:
                throw std::runtime_error(std::string(
                        "rotation can only take values in {0, 90, 180, 270}"
                ));
        }
    }

    Eigen::VectorXd Bicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (bicop_->rotation_) {
            case 0:
                return bicop_->hinv1(bicop_->cut_and_rotate(u));

            case 90:
                return bicop_->hinv2(bicop_->cut_and_rotate(u));

            case 180:
                return 1.0 - bicop_->hinv1(bicop_->cut_and_rotate(u)).array();

            case 270:
                return 1.0 - bicop_->hinv2(bicop_->cut_and_rotate(u)).array();

            default:
                throw std::runtime_error(std::string(
                        "rotation only takes value in {0, 90, 180, 270}"
                ));
        }
    }

    Eigen::VectorXd Bicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (bicop_->rotation_) {
            case 0:
                return bicop_->hinv2(bicop_->cut_and_rotate(u));

            case 90:
                return 1.0 - bicop_->hinv1(bicop_->cut_and_rotate(u)).array();

            case 180:
                return 1.0 - bicop_->hinv2(bicop_->cut_and_rotate(u)).array();

            case 270:
                return bicop_->hinv1(bicop_->cut_and_rotate(u));

            default:
                throw std::runtime_error(std::string(
                        "rotation only takes value in {0, 90, 180, 270}"
                ));
        }
    }
    //! @}


    //! Simulate from a bivariate copula
    //!
    //! @param n number of observations.
    //! @return Samples from the copula model.
    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::simulate(const int& n)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> U =
                tools_stats::simulate_uniform(n, 2);
        // use inverse Rosenblatt transform to generate a sample from the copula
        U.col(1) = hinv1(U);
        return U;
    }

    //! Fit statistics
    //!
    //! @param u \f$m \times 2\f$ matrix of observations.
    //! @{
    double Bicop::loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        Eigen::VectorXd ll = pdf(u);
        ll = ll.array().log();
        return ll.sum();
    }

    double Bicop::aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * loglik(u) + 2 * calculate_npars();
    }

    double Bicop::bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * loglik(u) + calculate_npars() * log(u.rows());
    }
    //! @}

    Eigen::MatrixXd Bicop::tau_to_parameters(const double& tau)
    {
        return bicop_->tau_to_parameters(tau);
    }

    double Bicop::parameters_to_tau(const Eigen::VectorXd& parameters)
    {
        return bicop_->parameters_to_tau(parameters);
    }


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
        return bicop_->get_rotation();
    }
    
    Eigen::MatrixXd Bicop::get_parameters() const 
    {
        return bicop_->get_parameters();
    }

    void Bicop::set_rotation(const int& rotation) {
        bicop_->set_rotation(rotation);
    }

    void Bicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        bicop_->set_parameters(parameters);
    }

    double Bicop::calculate_npars()
    {
        return bicop_->calculate_npars();
    }

    void Bicop::flip()
    {
        bicop_->flip();
    }

    BicopPtr Bicop::get_bicop()
    {
        return bicop_;
    };

    void Bicop::fit(
            const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
            std::string method
    )
    {
        bicop_->fit(data,method);
    }
    void Bicop::select(
            Eigen::Matrix<double, Eigen::Dynamic, 2> data,
            std::vector<BicopFamily> family_set,
            std::string method,
            std::string selection_criterion,
            bool preselect_families
    )
    {
        Bicop bicop(data, family_set, method,
                    selection_criterion, preselect_families);
        bicop_ = bicop.get_bicop();
    }
}
