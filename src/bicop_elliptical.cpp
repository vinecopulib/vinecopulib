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

#include "bicop_elliptical.hpp"


VecXd EllipticalBicop::hfunc2_default(const MatXd& u)
{
    return hfunc1_default(swap_cols(u));
}

VecXd EllipticalBicop::hinv2_default(const MatXd& u)
{
    return hinv1_default(swap_cols(u));
}

double EllipticalBicop::parameters_to_tau(const VecXd& parameters)
{
    double tau = (2 / M_PI) * asin(parameters(0));
    return tau;
}

// link between Kendall's tau and the par_bicop parameter
VecXd EllipticalBicop::tau_to_parameters(const double& tau)
{
    VecXd parameters = this->parameters_;
    parameters(0) = sin(tau * M_PI / 2);
    return parameters;
}

void EllipticalBicop::flip()
{
    // nothing to do because elliptical copulas are radially syemmetric
}
