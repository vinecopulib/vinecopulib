// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <cmath>

//! Shared helpers for the unit tests.
namespace test_utils {

//! @brief Element-wise comparison of two vectors that tolerates non-finite
//! values.
//!
//! Finite entries must agree within `atol + rtol * |b(i)|`, while non-finite
//! entries must agree exactly (both NaN, or equal infinities). This is useful
//! because some copula quantities are genuinely NaN/Inf at specific parameter
//! values (e.g., the Frank density at `theta = 0`, which is `0/0` and evaluates
//! to NaN or a finite value depending on the platform's floating-point
//! contraction), whereas `Eigen::DenseBase::isApprox` treats any NaN as a
//! mismatch. On failure the offending index and values are reported.
//!
//! @param a,b The vectors to compare.
//! @param rtol Relative tolerance for finite entries.
//! @param atol Absolute tolerance for finite entries.
inline ::testing::AssertionResult
all_close(const Eigen::VectorXd& a,
          const Eigen::VectorXd& b,
          double rtol = 1e-8,
          double atol = 1e-8)
{
  if (a.size() != b.size()) {
    return ::testing::AssertionFailure()
           << "size mismatch: " << a.size() << " != " << b.size();
  }
  for (Eigen::Index i = 0; i < a.size(); ++i) {
    const double x = a(i);
    const double y = b(i);
    if ((std::isfinite)(x) && (std::isfinite)(y)) {
      const double bound = atol + rtol * std::abs(y);
      if (std::abs(x - y) > bound) {
        return ::testing::AssertionFailure()
               << "mismatch at index " << i << ": " << x << " vs " << y
               << ", abs diff = " << std::abs(x - y)
               << ", tolerance = " << bound;
      }
    } else {
      // both NaN, or equal infinities (handles +inf/+inf and -inf/-inf)
      const bool both_nan = (std::isnan)(x) && (std::isnan)(y);
      if (!both_nan && !(x == y)) {
        return ::testing::AssertionFailure() << "non-finite mismatch at index "
                                             << i << ": " << x << " vs " << y;
      }
    }
  }
  return ::testing::AssertionSuccess();
}

} // namespace test_utils
