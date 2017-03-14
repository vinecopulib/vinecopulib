// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_class.hpp"
#include "bicop.hpp"
#include <exception>
#include <cmath>
#include "tools_stl.hpp"
#include "tools_stats.hpp"


//! Create a bivariate copula using the default contructor
//!
//! @param family the copula family.
//! @param rotation the rotation type.
//! @return A pointer to an object that inherits from \c Bicop.
BicopPtr Bicop::create(const int& family, const int& rotation)
{
    BicopPtr my_bicop;
    switch (family) {
        case 0:
            my_bicop = BicopPtr(new IndepBicop());
            break;
        case 1:
            my_bicop = BicopPtr(new GaussBicop());
            break;
        case 2:
            my_bicop = BicopPtr(new StudentBicop());
            break;
        case 3:
            my_bicop = BicopPtr(new ClaytonBicop());
            break;
        case 4:
            my_bicop = BicopPtr(new GumbelBicop());
            break;
        case 5:
            my_bicop = BicopPtr(new FrankBicop());
            break;
        case 6:
            my_bicop = BicopPtr(new JoeBicop());
            break;
        case 7:
            my_bicop = BicopPtr(new Bb1Bicop());
            break;
        case 8:
            my_bicop = BicopPtr(new Bb6Bicop());
            break;
        case 9:
            my_bicop = BicopPtr(new Bb7Bicop());
            break;
        case 10:
            my_bicop = BicopPtr(new Bb8Bicop());
            break;
        case 1001:
            my_bicop =  BicopPtr(new TrafokernelBicop());
            break;

        default:
            throw std::runtime_error(std::string("Family not implemented"));
    }

    my_bicop->set_rotation(rotation);
    return my_bicop;
}

//! Create a bivariate copula with a specified parameters vector
//!
//! @param family the copula family.
//! @param par the copula parameters (must be compatible with family).
//! @param rotation the rotation type.
//! @return A pointer to an object that inherits from \c Bicop.
BicopPtr Bicop::create(const int& family, const VecXd& parameters, const int& rotation)
{
    BicopPtr my_bicop = create(family, rotation);
    my_bicop->set_parameters(parameters);
    return my_bicop;
}

//! Select a bivariate copula
//!
//! @param data the data to fit the bivariate copula.
//! @param selection_criterion the selection criterion ("aic" or "bic").
//! @param family_set the set of copula families to consider (if empty, then 
//!     all families are included).
//! @param use_rotations whether rotations in the familyset are included.
//! @param preselect_families whether to exclude families before fitting based 
//!     on symmetry properties of the data.
//! @param method indicates the estimation method: either maximum likelihood 
//!     estimation (method = "mle") or inversion of Kendall's tau (method = 
//!     "itau"). When method = "itau" is used with families having more than one
//!      parameter, the main dependence parameter is found by inverting the 
//!      Kendall's tau and the remainders by a profile likelihood optimization.
//! @return A pointer to an object that inherits from \c Bicop.
BicopPtr Bicop::select(const MatXd& data,
                     std::string selection_criterion,
                     std::vector<int> family_set,
                     bool use_rotations,
                     bool preselect_families,
                     std::string method)
{
    std::vector<int> all_families = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1001};
    std::vector<int> itau_families = {0, 1, 2, 3, 4, 5, 6};
    std::vector<int> rotationless_families = {0, 1, 2, 5, 1001};

    // If the familyset is empty, use all families.
    // If the familyset is not empty, check that all included families are implemented.
    if (family_set.empty())
    {
        if (method == "itau")
        {
            family_set = tools_stl::intersect(all_families, itau_families);
        }
        else
        {
            family_set = all_families;
        }
    } else
    {
        for (unsigned int j = 0; j < family_set.size(); j++)
        {
            if (!tools_stl::is_member(family_set[j], all_families))
                throw std::runtime_error(std::string("One of the families is not implemented"));
        }
        if (method == "itau") {
            family_set = tools_stl::intersect(family_set, itau_families);
            if (family_set.empty())
                throw std::runtime_error(std::string("None of the families has method itau available"));
        }

    }

    auto temp_data = data;
    auto tau = pairwise_ktau(temp_data);

    // When using rotations, add only the ones that yield the appropriate 
    // association direction.
    std::vector<int> which_rotations = {0};
    if (use_rotations)
    {
        if (tau < 0)
        {
            which_rotations.pop_back();
            which_rotations.push_back(90);
            which_rotations.push_back(270);
        } else {
            which_rotations.push_back(180);
        }
    }

    double c1 = 0;
    double c2 = 0;
    if (preselect_families)
    {
        std::vector<double> c = get_c1c2(temp_data, tau);
        c1 = c[0];
        c2 = c[1];
    }

    // Create the combinations of families and rotations to estimate
    std::vector<int> families;
    std::vector<int> rotations;
    for (unsigned int j = 0; j < family_set.size(); j++)
    {
        bool is_rotationless = tools_stl::is_member(family_set[j], rotationless_families);
        bool preselect = true;
        if (is_rotationless)
        {
            if (preselect_families)
            {
                preselect = preselect_family(c1, c2, tau, family_set[j], 0, is_rotationless);
            }

            if (preselect)
            {
                families.push_back(family_set[j]);
                rotations.push_back(0);
            }
        } else {
            for (unsigned int k = 0; k < which_rotations.size(); k++)
            {
                if (preselect_families)
                {
                    preselect = preselect_family(c1, c2, tau, family_set[j], which_rotations[k], is_rotationless);
                }

                if (preselect)
                {
                    families.push_back(family_set[j]);
                    rotations.push_back(which_rotations[k]);
                }
            }
        }
    }

    // Estimate all models and select the best one using the selection_criterion
    BicopPtr fitted_bicop;
    double fitted_criterion = 1e6;
    for (unsigned int j = 0; j < families.size(); j++)
    {
        // Estimate the model
        BicopPtr new_bicop = create(families[j], rotations[j]);
        new_bicop->fit(temp_data, method);

        // Compute the selection criterion
        double new_criterion;
        if (selection_criterion == "aic")
        {
            new_criterion = new_bicop->aic(temp_data);
        }
        else if (selection_criterion == "bic")
        {
            new_criterion = new_bicop->bic(temp_data);
        }
        else
        {
            throw std::runtime_error(std::string("Selection criterion not implemented"));
        }

        // If the new model is better than the current one, then replace the current model by the new one
        if (new_criterion < fitted_criterion)
        {
            fitted_criterion = new_criterion;
            fitted_bicop = new_bicop;
        }
     }
    return fitted_bicop;

}

std::vector<double> get_c1c2(const MatXd& data, double tau)
{
    int n = data.rows();
    MatXd x = MatXd::Zero(n,2);
    MatXd z1 = x;
    MatXd z2 = x;
    x = qnorm(data);
    int count1 = 0, count2 = 0;
    for (int j = 0; j < n; ++j) {
        if (tau > 0)
        {
            if (x(j,0) > 0 && x(j,1) > 0)
            {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if (x(j,0) < 0 && x(j,1) < 0)
            {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        } else {
            if (x(j,0) < 0 && x(j,1) > 0)
            {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if (x(j,0) > 0 && x(j,1) < 0)
            {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        }
    }

    std::vector<double> c = {
        pairwise_cor(z1.block(0,0,count1-1,2)),
        pairwise_cor(z2.block(0,0,count2-1,2))
    };
    return c;
}

//! Copula density
//!
//! @param u \f$m \times 2\f$ matrix of evaluation points.
//! @return The copula density evaluated at \c u.
VecXd Bicop::pdf(const MatXd& u)
{
    VecXd f = pdf_default(cut_and_rotate(u));
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
VecXd Bicop::hfunc1(const MatXd& u)
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

VecXd Bicop::hfunc2(const MatXd& u)
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

VecXd Bicop::hinv1(const MatXd& u)
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

VecXd Bicop::hinv2(const MatXd& u)
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
MatXd Bicop::simulate(const int& n)
{
    MatXd U = simulate_uniform(n, 2);
    // use inverse Rosenblatt transform to generate a sample from the copula
    U.col(1) = hinv1(U);
    return U;
}

//! Fit statistics
//!
//! @param u \f$m \times 2\f$ matrix of observations.
//! @{
double Bicop::loglik(const MatXd& u)
{
    VecXd ll = this->pdf(u);
    ll = ll.array().log();
    return ll.sum();
}

double Bicop::aic(const MatXd& u)
{
    return -2 * this->loglik(u) + 2 * calculate_npars();
}

double Bicop::bic(const MatXd& u)
{
    return -2 * this->loglik(u) + this->calculate_npars() * log(u.rows());
}
//! @}


// TODO: generic fall-back that calculates Kendall's tau based on the pdf:
// tau = int_0^1 int_0^1 C(u, v) c(u, v) du dv
//     = int_0^1 int_0^1 (int_0^u int_0^v c(s, t) ds dt) c(u, v) du dv
double Bicop::calculate_tau()
{
    return 999.0;
}

VecXd Bicop::tau_to_parameters(const double& tau)
{
    throw std::runtime_error("Method not implemented for this family");
}


//! Getters and setters.
//! @{
int Bicop::get_family() const {return family_;}
int Bicop::get_rotation() const {return rotation_;}
std::string Bicop::get_association_direction() const {return association_direction_;}
VecXd Bicop::get_parameters() const {return parameters_;}
MatXd Bicop::get_parameters_bounds() const {return parameters_bounds_;}

void Bicop::set_rotation(const int& rotation) {
    check_rotation(rotation);
    rotation_ = rotation;
    if (this->association_direction_ == "positive" || this->association_direction_ == "negative")
    {
        if (rotation == 0 || rotation == 180)
        {
            this->association_direction_ = "positive";
        }
        else
        {
            this->association_direction_ = "negative";
        }
    }
}

void Bicop::set_parameters(const VecXd& parameters)
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
VecXd Bicop::hinv1_num(const MatXd &u)
{
    MatXd u_new = u;
    auto h1 = [&](const VecXd &v) {
        u_new.col(1) = v;
        return hfunc1_default(u_new);
    };
    return invert_f(u.col(1), h1);
}

VecXd Bicop::hinv2_num(const MatXd &u)
{
    MatXd u_new = u;
    auto h1 = [&](const VecXd &x) {
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
MatXd Bicop::cut_and_rotate(const MatXd& u)
{
    MatXd u_new(u.rows(), 2);

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
    MatXd eps = MatXd::Constant(u.rows(), 2, 1-10);
    u_new = u_new.array().min(1.0 - eps.array());
    u_new = u_new.array().max(eps.array());

    return u_new;
}


MatXd Bicop::swap_cols(const MatXd& u)
{
    MatXd u_swapped = u;
    u_swapped.col(0).swap(u_swapped.col(1));
    return u_swapped;
}
//! @}

//! Sanity checks
//! @{
void Bicop::check_parameters(const VecXd& parameters)
{
    int num_pars = parameters_bounds_.rows();

    if (parameters.size() != num_pars) {
        std::stringstream message;
        message <<
            "Wrong size of parameters for " << family_name_ << " copula; " <<
            "expected: " << num_pars << ", " <<
            "actual: " << parameters.size() << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    if (num_pars > 0) {
        std::stringstream message;
        for (int i = 0; i < num_pars; ++i) {
            if (parameters(i) < parameters_bounds_(i, 0)) {
                message <<
                    "parameters[" << i << "]" <<
                    " must be larger than " << parameters_bounds_(i, 0) <<
                    " for the " << family_name_ << " copula" << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
            if (parameters(i) > parameters_bounds_(i, 1)) {
                message <<
                "parameters[" << i << "]" <<
                    " must be smaller than " << parameters_bounds_(i, 1) <<
                    " for the " << family_name_ << " copula";
                throw std::runtime_error(message.str().c_str());
            }
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

bool preselect_family(double c1, double c2, double tau, int family, int rotation, bool is_rotationless)
{
    bool preselect = false;
    if (is_rotationless)
    {
        if (std::fabs(c1-c2) > 0.3)
        {
            if (family == 2 || family == 0 || family == 1001)
            {
                preselect = true;
            }
        } else
        {
            preselect = true;
        }
    } else
    {
        bool is_90or180 = (rotation == 90 || rotation == 180);
        if (family == 7 || family == 8 || family == 9 || family == 10)
        {
            if (tau > 0 && (rotation == 0 || rotation == 180))
            {
                preselect = true;
            }
            if (tau < 0 && (rotation == 90 || rotation == 270))
            {
                preselect = true;
            }
        }
        if (c1 - c2 > 0.05)
        {
            if (family == 3)
            {
                if (is_90or180)
                {
                    preselect = true;
                }
            }
            if (family == 4 || family == 6)
            {
                if (!is_90or180)
                {
                    preselect = true;
                }
            }
        } else if (c1 - c2 < -0.05)
        {
            if (family == 3)
            {
                if (!is_90or180)
                {
                    preselect = true;
                }
            }
            if (family == 4 || family == 6)
            {
                if (is_90or180)
                {
                    preselect = true;
                }
            }
        } else {
            if (tau > 0 && (rotation == 0 || rotation == 180))
            {
                preselect = true;
            }
            if (tau < 0 && (rotation == 90 || rotation == 270))
            {
                preselect = true;
            }
        }
    }
    return preselect;
}
