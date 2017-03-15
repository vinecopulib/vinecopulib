// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_indep.hpp"

namespace vinecopulib
{
    IndepBicop::IndepBicop()
    {
        family_ = 0;
        family_name_ = "Independence";
        rotation_ = 0;
        association_direction_ = "none";
    }

    IndepBicop::IndepBicop(const Eigen::VectorXd& parameters)
    {
        IndepBicop();
        set_parameters(parameters);
    }

    IndepBicop::IndepBicop(const Eigen::VectorXd& parameters, const int& rotation)
    {
        IndepBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

// PDF
    Eigen::VectorXd IndepBicop::pdf_default(const Eigen::MatrixXd& u)
    {
        return Eigen::VectorXd::Ones(u.rows());
    }

// hfunctions: the conditioning variable is put second
    Eigen::VectorXd IndepBicop::hfunc1_default(const Eigen::MatrixXd& u)
    {
        return u.col(1);
    }

    Eigen::VectorXd IndepBicop::hfunc2_default(const Eigen::MatrixXd& u)
    {
        return u.col(0);
    }

    Eigen::VectorXd IndepBicop::hinv1_default(const Eigen::MatrixXd& u)
    {
        return u.col(1);
    }

    Eigen::VectorXd IndepBicop::hinv2_default(const Eigen::MatrixXd& u)
    {
        return u.col(0);
    }

    Eigen::VectorXd IndepBicop::tau_to_parameters(const double &)
    {
        Eigen::VectorXd pars;
        return pars;
    }

    double IndepBicop::parameters_to_tau(const Eigen::VectorXd &)
    {
        return 0.0;
    }

    Eigen::VectorXd IndepBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }

    void IndepBicop::flip()
    {
        // nothing to do because independence copula is radially syemmetric
    }
}
