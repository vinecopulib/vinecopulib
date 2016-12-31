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

#include <exception>
#include "include/bicop.hpp"

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
            family_exists = (
                    std::find(
                            all_families.begin(),
                            all_families.end(),
                            family_set[j]
                    )  != all_families.end()
            );
            if (!family_exists)
                throw std::runtime_error(std::string("One of the families is not implemented"));
        }
    }

    // When using rotations, add only the ones that yield the appropriate association direction.
    std::vector<int> which_rotations = {0};
    if (use_rotations)
    {
        int n = data.rows();
        int d = 2;
        double tau = 0.0;
        MatXd newdata = data;
        ktau_matrix(newdata.data(), &d, &n, &tau);

        if (tau < 0)
        {
            which_rotations.pop_back();
            which_rotations.push_back(90);
            which_rotations.push_back(270);
        } else {
            which_rotations.push_back(180);
        }
        //std::cout << which_rotations[0] << "/" << which_rotations[1] << "/" << tau << std::endl;
    }

    // Create the combinations of families and rotations to estimate
    std::vector<int> families;
    std::vector<int> rotations;
    for (unsigned int j = 0; j < family_set.size(); j++)
    {
        bool is_rotationless = (
                std::find(
                        rotationless_families.begin(),
                        rotationless_families.end(),
                        family_set[j]
                )  != rotationless_families.end()
        );

        if (preselect_families)
        {
            // TODO
            /*MatXd x = MatXd::Zero(n,2);
            MatXd z = x;
            double c1, c2;
            for (int j = 0; j < n; ++j) {
                x(j,0) = gsl_cdf_ugaussian_Pinv(data(j, 0));
                x(j,1) = gsl_cdf_ugaussian_Pinv(data(j, 1));
            }
            if (tau > 0)
            {
                //std::cout << (x.array() > 0.0).rowwise().prod() << std::endl;

            }*/

        } else {
            if (is_rotationless)
            {
                families.push_back(family_set[j]);
                rotations.push_back(0);
            } else {
                for (unsigned int k = 0; k < which_rotations.size(); k++)
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
        //std::cout << families[j] << "/" << rotations[j] <<  std::endl;
        // Estimate the model
        BicopPtr new_bicop = create(families[j], rotations[j]);
        new_bicop->fit(data, method);

        // Compute the selection criterion
        double new_criterion;
        if (selection_criterion.compare("aic") == 0)
        {
            new_criterion = new_bicop->aic(data);
        } else if (selection_criterion.compare("bic") == 0)
        {
            new_criterion = new_bicop->bic(data);
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
        //std::cout << families[j] << " " << rotations[j] << " " << new_criterion << " " << new_bicop->loglik(data) << std::endl;
    }

    return fitted_bicop;

}

VecXd Bicop::pdf(const MatXd& u)
{
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
    std::default_random_engine generator(time(0));
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
    return -2 * this->loglik(u) + 2 * log(this->calculate_npars());
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
            fm(j) = hfunc1(v.row(j))(0) - u1(j);

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
            fm(j) = hfunc2(v.row(j))(0) - u1(j);

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