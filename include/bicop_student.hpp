// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_elliptical.hpp"

namespace vinecopulib
{
    class StudentBicop : public EllipticalBicop {

    public:
        // constructor
        StudentBicop();
        StudentBicop(const VectorXd& parameters);
        StudentBicop(const VectorXd& parameters, const int& rotation);

        // PDF
        VectorXd pdf_default(const MatrixXd& u);

        // hfunction
        VectorXd hfunc1_default(const MatrixXd& u);

        // inverse hfunction
        VectorXd hinv1_default(const MatrixXd& u);

    private:
        VectorXd get_start_parameters(const double tau);
    };
}
