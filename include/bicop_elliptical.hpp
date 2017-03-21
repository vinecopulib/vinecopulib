// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_parametric.hpp"

namespace vinecopulib
{
    class EllipticalBicop : public ParBicop
    {
    private:
        // hfunction and its inverse
        Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u);
        Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u);

        // link between Kendall's tau and the par_bicop parameter
        double parameters_to_tau(const Eigen::VectorXd& parameters);

        void flip();
    };
}
