// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "bicop_kernel.hpp"

namespace vinecopulib
{
    class TrafokernelBicop : public KernelBicop
    {
    public:
        TrafokernelBicop();
    private:
        void fit(const Eigen::MatrixXd& data, std::string);
    };
}
