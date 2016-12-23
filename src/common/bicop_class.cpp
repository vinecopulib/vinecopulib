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
#include "include/bicop_class.hpp"

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

double Bicop::loglik(MatXd& u)
{
    VecXd ll = this->pdf(u);
    ll = ll.array().log();
    return ll.sum();
}

double Bicop::aic(MatXd& u)
{
    double out = -2 * this->loglik(u) + 2 * calculate_npars();
    return out;
}

double Bicop::bic(MatXd& u)
{
    double out = -2 * this->loglik(u) + 2 * log(this->calculate_npars());
    return out;
}

// TODO: generic fall-back that calculates Kendall's tau based on the pdf:
// tau = int_0^1 int_0^1 C(u, v) c(u, v) du dv
//     = int_0^1 int_0^1 (int_0^u int_0^v c(s, t) ds dt) c(u, v) du dv
double Bicop::calculate_tau()
{
    return 999.0;
}
