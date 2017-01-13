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
    std::vector<int> all_families = {0, 1, 2, 3, 4, 5, 6};
    std::vector<int> rotationless_families = {0, 1, 2, 5};

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
    int count1 = 0, count2 = 0;
    for (int j = 0; j < n; ++j) {
        x(j,0) = gsl_cdf_ugaussian_Pinv(data(j, 0));
        x(j,1) = gsl_cdf_ugaussian_Pinv(data(j, 1));
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
    if (rotation_ == 0)
        return pdf_default(u);
    else
        return pdf_default(rotate_u(u));
}


VecXd Bicop::hfunc1(const MatXd& u)
{
    switch (rotation_) {
        case 0:
            return hfunc1_default(u);

        case 90:
            return hfunc2_default(rotate_u(u));

        case 180:
            return 1.0 - hfunc1_default(rotate_u(u)).array();

        case 270:
            return 1.0 - hfunc2_default(rotate_u(u)).array();

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
            return hfunc2_default(u);

        case 90:
            return 1.0 - hfunc1_default(rotate_u(u)).array();

        case 180:
            return 1.0 - hfunc2_default(rotate_u(u)).array();

        case 270:
            return hfunc1_default(rotate_u(u));

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
            return hinv1_default(u);

        case 90:
            return hinv2_default(rotate_u(u));

        case 180:
            return 1.0 - hinv1_default(rotate_u(u)).array();

        case 270:
            return 1.0 - hinv2_default(rotate_u(u)).array();

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
            return hinv2_default(u);

        case 90:
            return 1.0 - hinv1_default(rotate_u(u)).array();

        case 180:
            return 1.0 - hinv2_default(rotate_u(u)).array();

        case 270:
            return hinv1_default(rotate_u(u));

        default:
            throw std::runtime_error(std::string(
                "rotation only takes value in {0, 90, 180, 270}"
            ));
    }
}

MatXd Bicop::simulate(const int& n)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    MatXd U = MatXd::Zero(n, 2);
    for (int i = 0; i < n; ++i) {
        U(i, 0) = distribution(generator);
        U(i, 1) = distribution(generator);
    }
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

MatXd Bicop::rotate_u(const MatXd& u)
{
    MatXd u_rotated(u.rows(), u.cols());
    // counter-clockwise rotations
    switch (rotation_) {
        case 0:
            u_rotated = u;
            break;

        case 90:
            u_rotated.col(0) = u.col(1);
            u_rotated.col(1) = 1.0 - u.col(0).array();
            break;

        case 180:
            u_rotated.col(0) = 1.0 - u.col(0).array();
            u_rotated.col(1) = 1.0 - u.col(1).array();
            break;

        case 270:
            u_rotated.col(0) = 1.0 - u.col(1).array();
            u_rotated.col(1) = u.col(0);
            break;
    }

    return u_rotated;
}

MatXd Bicop::swap_cols(const MatXd& u)
{
    MatXd u_swapped = u;
    u_swapped.col(0).swap(u_swapped.col(1));
    return u_swapped;
}


VecXd Bicop::hinv1_num(const MatXd &u)
{
    MatXd v = u;
    VecXd u1 = u.col(1);
    int j = 0, br = 0, it = 0;
    double tol = 1e-12, xl = 1e-10, xh = 1-1e-10;

    v.col(1) = xl * VecXd::Ones(u1.size());
    VecXd fl = (hfunc1(v) - u1).cwiseAbs();
    v.col(1) = xh * VecXd::Ones(u1.size());
    VecXd fh = (hfunc1(v) - u1).cwiseAbs();
    VecXd fm = VecXd::Ones(u1.size());

    for (j = 0; j < u.rows(); ++j)
    {
        if (fl(j) <= tol) {
            v(j, 1) = xl;
            br = 1;
        }
        if (fh(j) <= tol) {
            v(j, 1) = xh;
            br = 1;
        }

        while (!br) {
            v(j, 1) = (xh + xl) / 2.0;
            fm(j) = hfunc1_default(v.row(j))(0) - u1(j);

            //stop if values become too close (avoid infinite loop)
            if (fabs(fm(j)) <= tol) br = 1;
            if (fabs(xl-xh) <= tol) br = 1;

            if (fm(j) < 0.0) {
                xl = v(j, 1);
                fh(j) = fm(j);
            } else {
                xh = v(j, 1);
                fl(j) = fm(j);
            }

            //stop if too many iterations are required (avoid infinite loop)
            ++it;
            if (it > 50) br = 1;
        }
        br = 0;
        it = 0;
        xl = 1e-10;
        xh = 1-1e-10;
    }

    return v.col(1);
}
VecXd Bicop::hinv2_num(const MatXd& u)
{
    MatXd v = u;
    VecXd u1 = u.col(0);
    int j = 0, br = 0, it = 0;
    double tol = 1e-12, xl = 1e-10, xh = 1-1e-10;

    v.col(0) = xl * VecXd::Ones(u1.size());
    VecXd fl = (hfunc2(v) - u1).cwiseAbs();
    v.col(0) = xh * VecXd::Ones(u1.size());
    VecXd fh = (hfunc2(v) - u1).cwiseAbs();
    VecXd fm = VecXd::Ones(u1.size());

    for (j = 0; j < u.rows(); ++j)
    {
        if (fl(j) <= tol) {
            v(j, 0) = xl;
            br = 1;
        }
        if (fh(j) <= tol) {
            v(j, 0) = xh;
            br = 1;
        }

        while (!br) {
            v(j, 0) = (xh + xl) / 2.0;
            fm(j) = hfunc2_default(v.row(j))(0) - u1(j);

            //stop if values become too close (avoid infinite loop)
            if (fabs(fm(j)) <= tol) br = 1;
            if (fabs(xl-xh) <= tol) br = 1;

            if (fm(j) < 0.0) {
                xl = v(j,0);
                fh(j) = fm(j);
            } else {
                xh = v(j,0);
                fl(j) = fm(j);
            }

            //stop if too many iterations are required (avoid infinite loop)
            ++it;
            if (it > 50) br = 1;
        }
        br = 0;
        it = 0;
        xl = 1e-10;
        xh = 1-1e-10;
    }

    return v.col(0);
}

void Bicop::set_rotation(const int& rotation) {
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
    std::stringstream message;

    if (parameters.size() != num_pars) {
        message << 
            "Wrong size of parameters for " << family_name_ << " copula; " <<
            "expected: " << num_pars << ", " <<
            "actual: " << parameters.size() << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    if (num_pars > 0) {
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

template<typename T> bool is_member(T element, std::vector<T> set)
{
    return std::find(set.begin(), set.end(), element) != set.end();
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
            if (family == 2 || family == 0)
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
