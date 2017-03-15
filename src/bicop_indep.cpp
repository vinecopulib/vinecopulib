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

    IndepBicop::IndepBicop(const VectorXd& parameters)
    {
        IndepBicop();
        set_parameters(parameters);
    }

    IndepBicop::IndepBicop(const VectorXd& parameters, const int& rotation)
    {
        IndepBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

// PDF
    VectorXd IndepBicop::pdf_default(const MatrixXd& u)
    {
        return VectorXd::Ones(u.rows());
    }

// hfunctions: the conditioning variable is put second
    VectorXd IndepBicop::hfunc1_default(const MatrixXd& u)
    {
        return u.col(1);
    }

    VectorXd IndepBicop::hfunc2_default(const MatrixXd& u)
    {
        return u.col(0);
    }

    VectorXd IndepBicop::hinv1_default(const MatrixXd& u)
    {
        return u.col(1);
    }

    VectorXd IndepBicop::hinv2_default(const MatrixXd& u)
    {
        return u.col(0);
    }

    VectorXd IndepBicop::tau_to_parameters(const double &)
    {
        VectorXd pars;
        return pars;
    }

    double IndepBicop::parameters_to_tau(const VectorXd &)
    {
        return 0.0;
    }

    VectorXd IndepBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }

    void IndepBicop::flip()
    {
        // nothing to do because independence copula is radially syemmetric
    }
}
