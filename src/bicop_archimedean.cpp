// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_archimedean.hpp"
#include "tools_stl.hpp"

namespace vinecopulib
{
    Eigen::VectorXd ArchimedeanBicop::pdf_default(const Eigen::MatrixXd& u)
    {
        auto f = [this](const double& u, const double& v) {
            double temp = generator_inv(generator(u) + generator(v));
            temp = generator_derivative2(temp)/std::pow(generator_derivative(temp), 3.0);
            return (-1)*generator_derivative(u)*generator_derivative(v)*temp;
        };
        return u.col(0).binaryExpr(u.col(1), f);
    }

    Eigen::VectorXd ArchimedeanBicop::hfunc1_default(const Eigen::MatrixXd& u)
    {
        auto f = [this](const double& u, const double& v) {
            double temp = generator_inv(generator(u) + generator(v));
            return generator_derivative(u)/generator_derivative(temp);
        };
        return u.col(0).binaryExpr(u.col(1), f);
    }

    Eigen::VectorXd ArchimedeanBicop::hfunc2_default(const Eigen::MatrixXd& u)
    {
        return hfunc1_default(swap_cols(u));
    }

    Eigen::VectorXd ArchimedeanBicop::hinv1_default(const Eigen::MatrixXd& u)
    {
        Eigen::VectorXd hinv = hinv1_num(u);
        return hinv;
    }

    Eigen::VectorXd ArchimedeanBicop::hinv2_default(const Eigen::MatrixXd& u)
    {
        return hinv1_default(swap_cols(u));
    }

    void ArchimedeanBicop::flip()
    {
        if (rotation_ == 90) {
            set_rotation(270);
        } else if (rotation_ == 270) {
            set_rotation(90);
        }
    }

    Eigen::VectorXd ArchimedeanBicop::get_start_parameters(const double)
    {
        Eigen::MatrixXd bounds = this->get_parameters_bounds();
        Eigen::VectorXd parameters = bounds.col(0) + Eigen::VectorXd::Constant(2, 0.1);
        return parameters;
    }

    double ArchimedeanBicop::flip_tau(double tau)
    {
        double flipped_tau = tau;
        if (tools_stl::is_member(rotation_, {90, 270})) {
            flipped_tau *= -1;
        }
        return flipped_tau;
    }
}
