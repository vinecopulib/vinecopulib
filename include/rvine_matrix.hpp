// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "tools_eigen.hpp"

namespace vinecopulib
{
    class RVineMatrix {
    public:
        //! \devgroup constructors Constructors
        //! @{
        RVineMatrix() {}
        RVineMatrix(const MatrixXi& matrix);
        //! @}

        //! Getters
        //! @{
        MatrixXi get_matrix();
        //! @}

        VectorXi get_order();
        MatrixXi in_natural_order();
        MatrixXi get_max_matrix();
        MatrixXb get_needed_hfunc1();
        MatrixXb get_needed_hfunc2();

        static MatrixXi construct_d_vine_matrix(const VectorXi& order);

    private:
        int d_;
        MatrixXi matrix_;
    };

    int relabel_one(const int& x, const VectorXi& old_labels, const VectorXi& new_labels);
    MatrixXi relabel_elements(const MatrixXi& matrix, const VectorXi& new_labels);
}
