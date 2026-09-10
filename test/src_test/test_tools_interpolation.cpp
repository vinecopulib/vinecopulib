// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_interpolation.hpp>

namespace test_tools_interpolation {
using namespace vinecopulib;
using tools_interpolation::InterpolationGrid;

namespace {

//! A grid whose margins are deliberately *not* exactly uniform, since it is
//! the per-grid-line rescaling of a non-uniform margin that the rectangle
//! probabilities have to reproduce.
InterpolationGrid
skewed_grid(int m)
{
  Eigen::VectorXd g(m);
  for (int i = 0; i < m; ++i) {
    g(i) = static_cast<double>(i) / (m - 1);
  }
  // a smooth positive surface, unequal spacing in the values
  Eigen::MatrixXd v(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      v(i, j) = 0.5 + std::exp(-3.0 * std::abs(g(i) - g(j))) + 0.2 * g(i);
    }
  }
  return InterpolationGrid(g, v, 0);
}

} // namespace

// Summing the rectangle probabilities over a partition of the first argument
// must leave the second argument's marginal increment, whatever the partition:
// the per-grid-line rescaling is the same in every term and the masses
// telescope. This is what makes the quantity a probability and not merely a
// mass, and it holds to the last few bits only because every weight is
// nonnegative -- a route through cumulative differences loses it.
TEST(tools_interpolation, rect_mass_telescopes_in_the_first_argument)
{
  auto grid = skewed_grid(30);
  for (int k : { 1, 2, 5, 16, 128, 1000 }) {
    for (double b2 : { 0.125, 0.5, 0.875, 1.0 }) {
      for (double width : { 1.0 / 8, 1.0 / 512 }) {
        const double a2 = std::max(b2 - width, 0.0);
        double total = 0.0;
        for (int i = 0; i < k; ++i) {
          total += grid.rect_mass(
            static_cast<double>(i) / k, static_cast<double>(i + 1) / k, a2, b2);
        }
        EXPECT_NEAR(total, b2 - a2, 1e-14)
          << "k = " << k << ", (a2, b2) = (" << a2 << ", " << b2 << ")";
      }
    }
  }
}

// And summing over a partition of the second argument that starts at zero must
// give the same answer as the single rectangle spanning it.
TEST(tools_interpolation, rect_mass_telescopes_in_the_second_argument)
{
  auto grid = skewed_grid(30);
  for (double x0 : { 0.0, 0.25 }) {
    for (double x1 : { 0.3, 1.0 }) {
      const double y = 0.75;
      const double whole = grid.rect_mass(x0, x1, 0.0, y);
      for (int k : { 3, 17, 256 }) {
        double total = 0.0;
        for (int i = 0; i < k; ++i) {
          total += grid.rect_mass(x0, x1, y * i / k, y * (i + 1) / k);
        }
        EXPECT_NEAR(total, whole, 1e-14)
          << "k = " << k << ", (x0, x1) = (" << x0 << ", " << x1 << ")";
      }
    }
  }
}

// A conditional distribution reaches 1, so a partition of the free argument
// sums to exactly 1 -- both the numerator and the denominator are sums of
// nonnegative terms, and the numerators partition the denominator.
TEST(tools_interpolation, cond_interval_mass_partitions_to_one)
{
  auto grid = skewed_grid(30);
  for (size_t cond_var : { 1u, 2u }) {
    for (double u_cond : { 0.02, 0.3, 0.5, 0.97 }) {
      for (int k : { 1, 4, 33, 512 }) {
        double total = 0.0;
        for (int i = 0; i < k; ++i) {
          total += grid.cond_interval_mass(u_cond,
                                           static_cast<double>(i) / k,
                                           static_cast<double>(i + 1) / k,
                                           cond_var);
        }
        EXPECT_NEAR(total, 1.0, 1e-14)
          << "cond_var = " << cond_var << ", u_cond = " << u_cond
          << ", k = " << k;
      }
    }
  }
}

// An empty or inverted interval carries no probability, and the bounds may
// arrive in either order -- rotating the data swaps a left limit past its own
// value, so the routines orient the rectangle themselves.
TEST(tools_interpolation, rect_mass_orients_its_own_bounds)
{
  auto grid = skewed_grid(30);
  const double p = grid.rect_mass(0.2, 0.35, 0.6, 0.9);
  EXPECT_GT(p, 0.0);
  EXPECT_DOUBLE_EQ(grid.rect_mass(0.35, 0.2, 0.6, 0.9), p);
  EXPECT_DOUBLE_EQ(grid.rect_mass(0.2, 0.35, 0.9, 0.6), p);
  EXPECT_DOUBLE_EQ(grid.rect_mass(0.35, 0.2, 0.9, 0.6), p);
  EXPECT_DOUBLE_EQ(grid.rect_mass(0.3, 0.3, 0.6, 0.9), 0.0);
  EXPECT_DOUBLE_EQ(grid.rect_mass(0.2, 0.35, 0.7, 0.7), 0.0);

  const double q = grid.cond_interval_mass(0.4, 0.2, 0.35, 1);
  EXPECT_GT(q, 0.0);
  EXPECT_DOUBLE_EQ(grid.cond_interval_mass(0.4, 0.35, 0.2, 1), q);
  EXPECT_DOUBLE_EQ(grid.cond_interval_mass(0.4, 0.3, 0.3, 1), 0.0);
}

// The exact route is the value the four-corner difference defines, so on a
// rectangle wide enough for that difference to be accurate the two agree.
TEST(tools_interpolation, rect_mass_agrees_with_a_cdf_difference)
{
  auto grid = skewed_grid(30);
  for (double a1 : { 0.1, 0.4 }) {
    for (double a2 : { 0.15, 0.55 }) {
      const double b1 = a1 + 0.3, b2 = a2 + 0.3;
      Eigen::MatrixXd corners(4, 2);
      corners << b1, b2, a1, b2, b1, a2, a1, a2;
      const Eigen::VectorXd c = grid.integrate_2d(corners);
      const double differenced = (c(0) + c(3)) - (c(1) + c(2));
      EXPECT_NEAR(grid.rect_mass(a1, b1, a2, b2), differenced, 1e-12)
        << "(a1, a2) = (" << a1 << ", " << a2 << ")";
    }
  }
}
}
