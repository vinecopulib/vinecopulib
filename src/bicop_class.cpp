/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "bicop.hpp"

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
        case 1001:
            my_bicop =  BicopPtr(new TrafokernelBicop());
            break;
        
        default:
            throw std::runtime_error(std::string("Family not implemented"));
    }
    
    my_bicop->set_rotation(rotation);
    return my_bicop;
}


BicopPtr Bicop::create(const int& family, const VecXd& parameters, const int& rotation)
{
    BicopPtr my_bicop = create(family, rotation);
    my_bicop->set_parameters(parameters);
    return my_bicop;
}

BicopPtr Bicop::select(const MatXd& data,
                     std::string selection_criterion,
                     std::vector<int> family_set,
                     bool use_rotations,
                     bool preselect_families,
                     std::string method)
{
    std::vector<int> all_families = {0, 1, 2, 3, 4, 5, 6, 1001};
    std::vector<int> rotationless_families = {0, 1, 2, 5, 1001};

    // If the familyset is empty, use all families.
    // If the familyset is not empty, check that all included families are implemented.
    if (family_set.empty())
    {
        family_set = all_families;
    } else
    {
        bool family_exists = true;
        for (unsigned int j = 0; j < family_set.size(); j++)
        {
            family_exists = is_member(family_set[j], all_families);
            if (!family_exists)
                throw std::runtime_error(std::string("One of the families is not implemented"));
        }
    }

    int n = data.rows();
    int d = 2;
    double tau = 0.0;
    MatXd newdata = data;
    ktau_matrix(newdata.data(), &d, &n, &tau);

    // When using rotations, add only the ones that yield the appropriate association direction.
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
        std::vector<double> c = get_c1c2(newdata, tau);
        c1 = c[0];
        c2 = c[1];
    }

    // Create the combinations of families and rotations to estimate
    std::vector<int> families;
    std::vector<int> rotations;
    for (unsigned int j = 0; j < family_set.size(); j++)
    {
        bool is_rotationless = is_member(family_set[j], rotationless_families);
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
        new_bicop->fit(newdata, method);
    
        // Compute the selection criterion
        double new_criterion;
        if (selection_criterion.compare("aic") == 0)
        {
            new_criterion = new_bicop->aic(newdata);
        } else if (selection_criterion.compare("bic") == 0)
        {
            new_criterion = new_bicop->bic(newdata);
        } else
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

    std::vector<double> c = {correlation(z1.block(0,0,count1-1,2)),correlation(z2.block(0,0,count2-1,2))};
    return c;
}

VecXd Bicop::pdf(const MatXd& u)
{
    return pdf_default(cut_and_rotate(u));
}


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

MatXd Bicop::simulate(const int& n)
{
    MatXd U = simulate_uniform(n, 2);
    // use inverse Rosenblatt transform to generate a sample from the copula
    U.col(1) = hinv1(U);
    return U;
}

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

// TODO: generic fall-back that calculates Kendall's tau based on the pdf:
// tau = int_0^1 int_0^1 C(u, v) c(u, v) du dv
//     = int_0^1 int_0^1 (int_0^u int_0^v c(s, t) ds dt) c(u, v) du dv
double Bicop::calculate_tau()
{
    return 999.0;
}

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
    MatXd eps = MatXd::Constant(u.rows(), 2, 1e-20);
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

void Bicop::set_rotation(const int& rotation) {
    check_rotation(rotation);    
    rotation_ = rotation;
    if ((this->association_direction_).compare("positive") == 0 ||
            (this->association_direction_).compare("negative") == 0)
    {
        if (rotation == 0 || rotation == 180)
        {
            this->association_direction_ = "positive";
        } else
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
    if (!is_member(rotation, allowed_rotations)) {
        std::string message = "rotation must be one of {0, 90, 180, 270}";
        throw std::runtime_error(message);
    }
        
}

double correlation(const MatXd& z)
{
    double rho;
    MatXd x = z.rowwise() - z.colwise().mean();
    MatXd sigma = x.adjoint() * x;
    rho = sigma(1,0) / sqrt(sigma(0,0) * sigma(1,1));

    return rho;
}


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

//! Numerical inversion of a function on [0, 1].
//! 
//! Computes the inverse \f$f^{-1}\f$ of a function \f$f\f$ by the bisection 
//! method.
//! 
//! @param x evaluation points.
//! @param f the function to invert.
//! @param n_iter the number of iterations for the bisection (defaults to 35,
//! guaranteeing an accuracy of 0.5^35 ~= 6e-11). 
//! 
//! @return f^{-1}(x).
VecXd invert_f(const VecXd &x, std::function<VecXd(const VecXd&)> f, int n_iter) 
{
    VecXd xl = VecXd::Constant(x.size(), 1e-20);
    VecXd xh = 1.0 - xl.array();
    VecXd x_tmp = x;
    for (int iter = 0; iter < n_iter; ++iter) {
        x_tmp = (xh + xl) / 2.0;
        VecXd fm = f(x_tmp) - x;
        xl = (fm.array() < 0).select(x_tmp, xl);
        xh = (fm.array() < 0).select(xh, x_tmp);
    }

    return x_tmp;
}
