/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "bicop_archimedean.hpp"

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
    VecXd hinv = hinv1_num(u);
    return hinv;
}

VecXd ArchimedeanBicop::hinv2_default(const MatXd& u)
{
    return hinv1_default(swap_cols(u));
}

void ArchimedeanBicop::flip()
{
    if (rotation_ == 90) {
        set_rotation(270);
    } else if (rotation_ == 270) {
        set_rotation(90);
    }
}

VecXd ArchimedeanBicop::get_start_parameters(const double tau)
{
    MatXd bounds = this->get_parameters_bounds();
    VecXd parameters = bounds.col(0) + VecXd::Constant(2, 0.1);
    return parameters;
}
