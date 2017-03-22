// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/class.hpp"

namespace vinecopulib
{
    Bicop::Bicop()
    {
        bicop_ = AbstractBicop::create();
    }

    Bicop::Bicop(BicopFamily family, int rotation)
    {
        bicop_ = AbstractBicop::create(family, rotation);
    }
    
    Bicop::Bicop(BicopFamily family, int rotation, Eigen::VectorXd parameters)
    {
        bicop_ = AbstractBicop::create(family, rotation, parameters);
    }

    Bicop::Bicop(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                 std::vector<BicopFamily> family_set,
                 std::string method,
                 std::string selection_criterion,
                 bool preselect_families)
    {
        bicop_ = AbstractBicop::select(data,
                                       family_set,
                                       method,
                                       selection_criterion,
                                       preselect_families);
    }

    Eigen::VectorXd Bicop::pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->pdf(u);
    }

    Eigen::VectorXd Bicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->hfunc1(u);
    }

    Eigen::VectorXd Bicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->hfunc2(u);
    }

    Eigen::VectorXd Bicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->hinv1(u);
    }

    Eigen::VectorXd Bicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->hinv2(u);
    }

    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::simulate(const int& n)
    {
        return bicop_->simulate(n);
    }

    double Bicop::loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->loglik(u);
    }

    double Bicop::aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->aic(u);
    }

    double Bicop::bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return bicop_->bic(u);
    }

    Eigen::MatrixXd Bicop::tau_to_parameters(const double& tau)
    {
        return bicop_->tau_to_parameters(tau);
    }


    BicopFamily Bicop::get_family() const 
    {
        return bicop_->get_family();
    }
    
    std::string Bicop::get_family_name() const 
    {
        return bicop_->get_family_name();
    };
    
    int Bicop::get_rotation() const
    {
        return bicop_->get_rotation();
    }
    
    Eigen::MatrixXd Bicop::get_parameters() const 
    {
        return bicop_->get_parameters();
    }

    void Bicop::set_rotation(const int& rotation) {
        bicop_->set_rotation(rotation);
    }

    void Bicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        bicop_->set_parameters(parameters);
    }

    void Bicop::flip()
    {
        bicop_->flip();
    }
}
