// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include <exception>
#include <cmath>

#include "bicop/bicop.hpp"
#include "bicop/class.hpp"
#include "bicop/tools_bicopselect.hpp"
#include "misc/tools_stl.hpp"
#include "misc/tools_stats.hpp"

namespace vinecopulib
{
    //! Create a bivariate copula using the default contructor
    //!
    //! @param family the copula family.
    //! @param rotation the rotation type.
    //! @param parameters the copula parameters (optional, must be compatible 
    //!     with family).
    //! @return A pointer to an object that inherits from \c Bicop.
    //! @{
    BicopPtr Bicop::create(BicopFamily family, int rotation)
    {
        BicopPtr new_bicop;
        switch (family) {
            case BicopFamily::indep:
                new_bicop = BicopPtr(new IndepBicop());
                break;
            case BicopFamily::gaussian:
                new_bicop = BicopPtr(new GaussianBicop());
                break;
            case BicopFamily::student:
                new_bicop = BicopPtr(new StudentBicop());
                break;
            case BicopFamily::clayton:
                new_bicop = BicopPtr(new ClaytonBicop());
                break;
            case BicopFamily::gumbel:
                new_bicop = BicopPtr(new GumbelBicop());
                break;
            case BicopFamily::frank:
                new_bicop = BicopPtr(new FrankBicop());
                break;
            case BicopFamily::joe:
                new_bicop = BicopPtr(new JoeBicop());
                break;
            case BicopFamily::bb1:
                new_bicop = BicopPtr(new Bb1Bicop());
                break;
            case BicopFamily::bb6:
                new_bicop = BicopPtr(new Bb6Bicop());
                break;
            case BicopFamily::bb7:
                new_bicop = BicopPtr(new Bb7Bicop());
                break;
            case BicopFamily::bb8:
                new_bicop = BicopPtr(new Bb8Bicop());
                break;
            case BicopFamily::tll0:
                new_bicop =  BicopPtr(new TrafokernelBicop());
                break;

            default:
                throw std::runtime_error(std::string("Family not implemented"));
        }
        
        new_bicop->set_rotation(rotation);
            
        return new_bicop;
    }
    
    BicopPtr Bicop::create(
        BicopFamily family,
        int rotation,
        Eigen::VectorXd parameters
    )
    {
        auto new_bicop = Bicop::create(family, rotation);
        new_bicop->set_parameters(parameters);
        return new_bicop;
    }
    
    //!@}

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
    BicopPtr Bicop::select(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
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
            bool is_rotationless = is_member(family, bicop_families::rotationless);
            bool preselect = true;
            if (is_rotationless) {
                if (preselect_families) {
                    preselect = preselect_family(c, tau, family, 0, is_rotationless);
                }
                if (preselect) {
                    families.push_back(family);
                    rotations.push_back(0);
                }
            } else {
                for (auto rotation : which_rotations) {
                    if (preselect_families) {
                        preselect = preselect_family(c, tau, family, rotation, is_rotationless);
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
            BicopPtr new_bicop = create(families[j], rotations[j]);
            new_bicop->fit(data, method);

            // Compute the selection criterion
            double new_criterion;
            if (selection_criterion == "aic") {
                new_criterion = new_bicop->aic(data);
            } else if (selection_criterion == "bic") {
                new_criterion = new_bicop->bic(data);
            } else {
                throw std::runtime_error(std::string("Selection criterion not implemented"));
            }

            // If the new model is better than the current one, then replace the current model by the new one
            if (new_criterion < fitted_criterion) {
                fitted_criterion = new_criterion;
                fitted_bicop = new_bicop;
            }
        }
        return fitted_bicop;

    }

    //! Copula density
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The copula density evaluated at \c u.
    Eigen::VectorXd Bicop::pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        Eigen::VectorXd f = pdf_default(cut_and_rotate(u));
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
    //! rotations are properly handled.  They call \c hfunc1_default,
    //! \c hfunc2_default, \c hinv1_default, and hinv2_default which are
    //! family-specific implementations for `rotation = 0`.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The (inverse) h-function evaluated at \c u.
    //! @{
    Eigen::VectorXd Bicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        switch (rotation_) {
            case 0:
                return hfunc1_default(cut_and_rotate(u));

            case 90:
                return hfunc2_default(cut_and_rotate(u));

            case 180:
                return 1.0 - hfunc1_default(cut_and_rotate(u)).array();

            case 270:
                return 1.0 - hfunc2_default(cut_and_rotate(u)).array();

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
                return hfunc2_default(cut_and_rotate(u));

            case 90:
                return 1.0 - hfunc1_default(cut_and_rotate(u)).array();

            case 180:
                return 1.0 - hfunc2_default(cut_and_rotate(u)).array();

            case 270:
                return hfunc1_default(cut_and_rotate(u));

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
                return hinv1_default(cut_and_rotate(u));

            case 90:
                return hinv2_default(cut_and_rotate(u));

            case 180:
                return 1.0 - hinv1_default(cut_and_rotate(u)).array();

            case 270:
                return 1.0 - hinv2_default(cut_and_rotate(u)).array();

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
                return hinv2_default(cut_and_rotate(u));

            case 90:
                return 1.0 - hinv1_default(cut_and_rotate(u)).array();

            case 180:
                return 1.0 - hinv2_default(cut_and_rotate(u)).array();

            case 270:
                return hinv1_default(cut_and_rotate(u));

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
        Eigen::VectorXd ll = this->pdf(u);
        ll = ll.array().log();
        return ll.sum();
    }

    double Bicop::aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * this->loglik(u) + 2 * calculate_npars();
    }

    double Bicop::bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * this->loglik(u) + this->calculate_npars() * log(u.rows());
    }
    //! @}

    Eigen::MatrixXd Bicop::tau_to_parameters(const double& tau)
    {
        return tau_to_parameters_default(tau);
    }

    Eigen::VectorXd no_tau_to_parameters(const double&)
    {
        throw std::runtime_error("Method not implemented for this family");
    }

    //! Getters and setters.
    //! @{
    BicopFamily Bicop::get_family() const 
    {
        return family_;
    }
    
    std::string Bicop::get_family_name() const 
    {
        return vinecopulib::get_family_name(family_);
    };
    
    int Bicop::get_rotation() const
    {
        return rotation_;
    }
    
    Eigen::MatrixXd Bicop::get_parameters() const 
    {
        return parameters_;
    }
    
    Eigen::MatrixXd Bicop::get_parameters_lower_bounds() const 
    {
        return parameters_lower_bounds_;
    }
    
    Eigen::MatrixXd Bicop::get_parameters_upper_bounds() const 
    {
        return parameters_upper_bounds_;
    }

    void Bicop::set_rotation(const int& rotation) {
        check_rotation(rotation);
        rotation_ = rotation;
    }

    void Bicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        check_parameters(parameters);
        parameters_ = parameters;
    }
    //! @}

    //! Numerical inversion of h-functions
    //!
    //! These are generic functions to invert the hfunctions numerically.
    //! They can be used in derived classes to define \c hinv1 and \c hinv2.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The numerical inverse of h-functions.
    //! @{
    Eigen::VectorXd Bicop::hinv1_num(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> u_new = u;
        auto h1 = [&](const Eigen::VectorXd &v) {
            u_new.col(1) = v;
            return hfunc1_default(u_new);
        };
        return invert_f(u.col(1), h1);
    }

    Eigen::VectorXd Bicop::hinv2_num(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> u_new = u;
        auto h1 = [&](const Eigen::VectorXd &x) {
            u_new.col(0) = x;
            return hfunc2_default(u_new);
        };

        return invert_f(u.col(0), h1);
    }
    //! @}

    //! Data manipulations for rotated families
    //!
    //! @param u \f$m \times 2\f$ matrix of data.
    //! @return The manipulated data.
    //! @{
    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::cut_and_rotate(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
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
            Eigen::Matrix<double, Eigen::Dynamic, 2>::Constant(u.rows(), 2, 1-10);
        u_new = u_new.array().min(1.0 - eps.array());
        u_new = u_new.array().max(eps.array());

        return u_new;
    }


    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::swap_cols(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> u_swapped = u;
        u_swapped.col(0).swap(u_swapped.col(1));
        return u_swapped;
    }
    //! @}

    
    //! Sanity checks
    //! @{
    void Bicop::check_parameters(const Eigen::MatrixXd& parameters)
    {
        check_parameters_size(parameters);
        check_parameters_lower(parameters);
        check_parameters_upper(parameters);
    }
    
    
    void Bicop::check_parameters_size(const Eigen::MatrixXd& parameters)
    {
        if (parameters.size() != parameters_.size()) {
            if (parameters.rows() != parameters_.rows()) {
                std::stringstream message;
                message <<
                    "parameters have has wrong number of rows " << 
                    "for " << get_family_name() << " copula; " << 
                    "expected: " << parameters_.rows() << ", " <<
                    "actual: " << parameters.rows() << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
            if (parameters.cols() != parameters_.cols()) {
                std::stringstream message;
                message <<
                    "parameters have wrong number of columns " << 
                    "for " << get_family_name() << " copula; " << 
                    "expected: " << parameters_.cols() << ", " <<
                    "actual: " << parameters.cols() << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
        }
    }
    
    
    void Bicop::check_parameters_lower(const Eigen::MatrixXd& parameters)
    {
        if (parameters_lower_bounds_.size() > 0) {
            std::stringstream message;
            if ((parameters.array() < parameters_lower_bounds_.array()).any()) {
                message <<
                    "parameters exceed lower bound " << 
                    " for " << get_family_name() << " copula; \n" << 
                    "bound: \n" << parameters_lower_bounds_ << "\n" <<
                    "actual: " << parameters << "\n";
                throw std::runtime_error(message.str().c_str());  
            }
        }
    }
    
    void Bicop::check_parameters_upper(const Eigen::MatrixXd& parameters)
    {
        if (parameters_upper_bounds_.size() > 0) {
            std::stringstream message;
            if ((parameters.array() > parameters_upper_bounds_.array()).any()) {
                message <<
                    "parameters exceed upper bound " << 
                    " for " << get_family_name() << " copula; \n" << 
                    "bound: \n" << parameters_upper_bounds_ << "\n" <<
                    "actual: " << parameters << "\n";
                throw std::runtime_error(message.str().c_str());  
            }
        }
    }

    void Bicop::check_rotation(const int& rotation)
    {
        std::vector<int> allowed_rotations = {0, 90, 180, 270};
        if (!tools_stl::is_member(rotation, allowed_rotations)) {
            std::string message = "rotation must be one of {0, 90, 180, 270}";
            throw std::runtime_error(message);
        }
    }
    //! @}
}
