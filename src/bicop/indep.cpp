// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/indep.hpp"

namespace vinecopulib
{
    IndepBicop::IndepBicop()
    {
        family_ = BicopFamily::indep;
        rotation_ = 0;
    }

    Eigen::VectorXd IndepBicop::pdf_default(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return Eigen::VectorXd::Ones(u.rows());
    }

    Eigen::VectorXd IndepBicop::hfunc1_default(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(1);
    }

    Eigen::VectorXd IndepBicop::hfunc2_default(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(0);
    }

    Eigen::VectorXd IndepBicop::hinv1_default(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(1);
    }

    Eigen::VectorXd IndepBicop::hinv2_default(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(0);
    }

    Eigen::MatrixXd IndepBicop::tau_to_parameters_default(const double &)
    {
        Eigen::VectorXd pars(0);
        return pars;
    }

    double IndepBicop::parameters_to_tau(const Eigen::VectorXd&)
    {
        return 0.0;
    }

    Eigen::VectorXd IndepBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters_default(tau);
    }

    void IndepBicop::flip()
    {
        // nothing to do because independence copula is radially syemmetric
    }
}
