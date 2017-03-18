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
        family_ = BicopFamily::indep;
        rotation_ = 0;
        parameters_ = {Eigen::VectorXd()};
        parameters_lower_bounds_ = {Eigen::VectorXd()};
        parameters_upper_bounds_ = {Eigen::VectorXd()};
    }

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

    std::vector<Eigen::MatrixXd> IndepBicop::tau_to_parameters(const double &)
    {
        return {Eigen::VectorXd()};
    }

    double IndepBicop::parameters_to_tau(const std::vector<Eigen::MatrixXd> &)
    {
        return 0.0;
    }

    std::vector<Eigen::MatrixXd> IndepBicop::get_start_parameters(const double tau)
    {
        return {tau_to_parameters(tau)};
    }

    void IndepBicop::flip()
    {
        // nothing to do because independence copula is radially syemmetric
    }
}
