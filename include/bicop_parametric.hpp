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

#ifndef VINECOPULIB_BICOP_PARAMETRIC_HPP
#define VINECOPULIB_BICOP_PARAMETRIC_HPP

#include "bicop_class.hpp"

class ParBicop : public Bicop {

public:
    // fit copula parameters
    void fit(const MatXd& data, std::string method);

    // link between Kendall's tau and the par_bicop parameter
    virtual double parameters_to_tau(const VecXd& parameters) = 0;
    double calculate_tau() {return this->parameters_to_tau(parameters_);}

    // number of parameters
    double calculate_npars();

private:
    virtual VecXd get_start_parameters(const double tau) = 0;
};

#endif
