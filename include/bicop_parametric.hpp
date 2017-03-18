// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_class.hpp"

namespace vinecopulib
{
    class ParBicop : public Bicop
    {
    public:
        // fit copula parameters
        void fit(const Eigen::MatrixXd& data, std::string method);

        // link between Kendall's tau and the par_bicop parameter
        virtual double parameters_to_tau(const std::vector<Eigen::MatrixXd>& parameters) = 0;
        double calculate_tau() {return this->parameters_to_tau(parameters_);}

        // number of parameters
        double calculate_npars();

    private:
        virtual std::vector<Eigen::MatrixXd> get_start_parameters(const double tau) = 0;
    };
}
