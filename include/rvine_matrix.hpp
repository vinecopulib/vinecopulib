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
        RVineMatrix(const MatXi& matrix);
        //! @}

        //! Getters
        //! @{
        MatXi get_matrix();
        //! @}

        VecXi get_order();
        MatXi in_natural_order();
        MatXi get_max_matrix();
        MatXb get_needed_hfunc1();
        MatXb get_needed_hfunc2();

        static MatXi construct_d_vine_matrix(const VecXi& order);

    private:
        int d_;
        MatXi matrix_;
    };

    int relabel_one(const int& x, const VecXi& old_labels, const VecXi& new_labels);
    MatXi relabel_elements(const MatXi& matrix, const VecXi& new_labels);
}
