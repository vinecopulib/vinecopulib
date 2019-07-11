// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <iostream>
#include <vinecopulib/misc/tools_integration_legendre.hpp>

namespace vinecopulib {

//! Utilities for integration
namespace tools_integration {

//! @brief load quadrature nodes and weights for numerical integration
//!
//! @param n number of nodes and weights (10, 20, 30, 40 or 50)
//!
//! @return An \f$ n \times 2 \f$ with the nodes in the first colum and the
//! weights in the second.
inline Eigen::Matrix<double, Eigen::Dynamic, 2>
quadrature_rule(const size_t n_nodes)
{
  Eigen::Matrix<double, Eigen::Dynamic, 2> output;
  switch (n_nodes) {
    case 10:
      output = tools_legendre::nodes_weights_10;
      break;
    case 20:
      output = tools_legendre::nodes_weights_20;
      break;
    case 30:
      output = tools_legendre::nodes_weights_30;
      break;
    case 40:
      output = tools_legendre::nodes_weights_40;
      break;
    case 50:
      output = tools_legendre::nodes_weights_50;
      break;
    default:
      throw std::runtime_error("n must be one of 10, 20, 30, 40, 50.");
  }
  return output;
}

//! @brief integrate a function.
//!
//! @param f the function to be integrated
//! @param lower lower integration bound
//! @param upper upper integration bound
//! @param n_nodes number of quadrature nodes and weights (10, 20, 30, 40 or 50)
//!
//! @return the value of the integral
inline double
integrate(std::function<double(double)> f,
          const double lower,
          const double upper,
          const size_t n_nodes)
{
  auto nodes_weights = quadrature_rule(n_nodes);
  nodes_weights.col(0) = lower + (upper - lower) * nodes_weights.col(0).array();
  auto f_vals = tools_eigen::unaryExpr_or_nan(nodes_weights.col(0), f);
  return (lower - upper) * nodes_weights.col(1).cwiseProduct(f_vals).sum();
}
}
}
