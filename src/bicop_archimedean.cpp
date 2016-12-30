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

#include "include/bicop_archimedean.hpp"


// PDF
VecXd ArchimedeanBicop::pdf_default(const MatXd& u)
{
    VecXd f = VecXd::Ones(u.rows());
    VecXd v = generator(u.col(0)) + generator(u.col(1));
    f = generator_derivative(u.col(0)).cwiseProduct(generator_derivative(u.col(1)));
    VecXd numerator = generator_derivative2(generator_inv(v));
    VecXd denominator = generator_derivative(generator_inv(v)).array().pow(3.0);
    f = (-1)*f.cwiseProduct(numerator);
    f = f.cwiseQuotient(denominator);

    return f;
}

VecXd ArchimedeanBicop::hfunc1_default(const MatXd& u)
{
    VecXd h(u.rows());
    VecXd v(u.rows());

    v = generator(u.col(0)) + generator(u.col(1));
    h = generator_derivative(u.col(0)).cwiseQuotient(generator_derivative(generator_inv(v)));

    return h;
}

VecXd ArchimedeanBicop::hfunc2_default(const MatXd& u)
{
    return hfunc1_default(swap_cols(u));
}

VecXd ArchimedeanBicop::hinv1_default(const MatXd& u)
{
    return hinv(u);
}

VecXd ArchimedeanBicop::hinv2_default(const MatXd& u)
{
    return hinv(swap_cols(u));
}
