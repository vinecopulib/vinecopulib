// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop/class.hpp"

namespace vinecopulib
{
    class ParBicop : public Bicop
    {
    private:
        void fit(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
             std::string method
         );
        double calculate_npars();
        virtual Eigen::VectorXd get_start_parameters(const double tau) = 0;
    };
}
