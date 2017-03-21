// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "tools_eigen.hpp"

namespace vinecopulib
{
    class RVineMatrix
    {
    public:
        //! \devgroup constructors Constructors
        //! @{
        RVineMatrix() {}
        RVineMatrix(const Eigen::MatrixXi& matrix);
        //! @}

        //! Getters
        //! @{
        Eigen::MatrixXi get_matrix() const;
        //! @}

        Eigen::VectorXi get_order() const;
        Eigen::MatrixXi in_natural_order() const;
        Eigen::MatrixXi get_max_matrix() const;
        MatrixXb get_needed_hfunc1() const;
        MatrixXb get_needed_hfunc2() const;

        static Eigen::MatrixXi construct_d_vine_matrix(const Eigen::VectorXi& order);

    private:
        int d_;
        Eigen::MatrixXi matrix_;
    };

    int relabel_one(const int& x, const Eigen::VectorXi& old_labels, const Eigen::VectorXi& new_labels);
    Eigen::MatrixXi relabel_elements(const Eigen::MatrixXi& matrix, const Eigen::VectorXi& new_labels);
}
