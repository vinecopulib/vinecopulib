// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <functional>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_integration {

Eigen::Matrix<double, Eigen::Dynamic, 2>
quadrature_rule(const size_t n_nodes = 20);

double
integrate(std::function<double(double)> f,
          const double lower = 0,
          const double upper = 1,
          const size_t n_nodes = 20);
}
}

#include <vinecopulib/misc/implementation/tools_integration.ipp>
