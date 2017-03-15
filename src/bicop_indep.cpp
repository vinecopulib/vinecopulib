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

    IndepBicop::IndepBicop(const VecXd& parameters)
    {
        IndepBicop();
        set_parameters(parameters);
    }

    IndepBicop::IndepBicop(const VecXd& parameters, const int& rotation)
    {
        IndepBicop();
        set_parameters(parameters);
        set_rotation(rotation);
    }

// PDF
    VecXd IndepBicop::pdf_default(const MatXd& u)
    {
        return VecXd::Ones(u.rows());
    }

// hfunctions: the conditioning variable is put second
    VecXd IndepBicop::hfunc1_default(const MatXd& u)
    {
        return u.col(1);
    }

    VecXd IndepBicop::hfunc2_default(const MatXd& u)
    {
        return u.col(0);
    }

    VecXd IndepBicop::hinv1_default(const MatXd& u)
    {
        return u.col(1);
    }

    VecXd IndepBicop::hinv2_default(const MatXd& u)
    {
        return u.col(0);
    }

    VecXd IndepBicop::tau_to_parameters(const double &)
    {
        VecXd pars;
        return pars;
    }

    double IndepBicop::parameters_to_tau(const VecXd &)
    {
        return 0.0;
    }

    VecXd IndepBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }

    void IndepBicop::flip()
    {
        // nothing to do because independence copula is radially syemmetric
    }
}
