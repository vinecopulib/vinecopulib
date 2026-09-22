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

// A narrow interval is built from the overlap's own width, not from a
// difference of cumulative integrals, so it stays exact whether it sits inside
// one cell or straddles a grid point -- and splitting it must reproduce it.
TEST(tools_interpolation, narrow_intervals_are_exact_within_and_across_cells)
{
  auto grid = skewed_grid(17); // grid points at k / 16
  const double knot = 8.0 / 16;
  for (size_t cond_var : { 1u, 2u }) {
    for (double w : { 1e-2, 1e-5, 1e-9 }) {
      for (double lo : { knot + 1e-2, knot - w / 2 }) { // inside, straddling
        const double m = grid.cond_interval_mass(0.4, lo, lo + w, cond_var);
        EXPECT_GT(m, 0.0);
        const double split =
          grid.cond_interval_mass(0.4, lo, lo + w / 2, cond_var) +
          grid.cond_interval_mass(0.4, lo + w / 2, lo + w, cond_var);
        EXPECT_NEAR(m, split, 1e-12 * m)
          << "cond_var " << cond_var << ", width " << w << ", lo " << lo;
      }
    }
  }
}

// An interval mass is not a clipped probability: the whole line is exactly one,
// an empty interval exactly zero, and bounds outside the unit interval clip to
// it rather than contributing.
TEST(tools_interpolation, conditional_interval_mass_handles_the_boundary)
{
  auto grid = skewed_grid(30);
  for (size_t cond_var : { 1u, 2u }) {
    EXPECT_NEAR(grid.cond_interval_mass(0.4, 0.0, 1.0, cond_var), 1.0, 1e-14);
    EXPECT_DOUBLE_EQ(grid.cond_interval_mass(0.4, 0.3, 0.3, cond_var), 0.0);
    EXPECT_DOUBLE_EQ(grid.cond_interval_mass(0.4, 1.0, 1.0, cond_var), 0.0);
    EXPECT_DOUBLE_EQ(grid.cond_interval_mass(0.4, -0.5, 0.25, cond_var),
                     grid.cond_interval_mass(0.4, 0.0, 0.25, cond_var));
    EXPECT_DOUBLE_EQ(grid.cond_interval_mass(0.4, 0.75, 1.5, cond_var),
                     grid.cond_interval_mass(0.4, 0.75, 1.0, cond_var));
  }
}

// `integrate_1d` is that mass over `[0, u]`, clipped because it is consumed as
// a probability, and `inverse_integrate_1d` inverts it.
TEST(tools_interpolation, conditional_cdf_and_quantile_agree)
{
  auto grid = skewed_grid(30);
  for (size_t cond_var : { 1u, 2u }) {
    for (double u_cond : { 0.1, 0.5, 0.9 }) {
      for (double u : { 0.05, 0.25, 0.5, 0.75, 0.95 }) {
        Eigen::MatrixXd q(1, 2);
        q << (cond_var == 1 ? u_cond : u), (cond_var == 1 ? u : u_cond);
        const double p = grid.integrate_1d(q, cond_var)(0);
        EXPECT_NEAR(p, grid.cond_interval_mass(u_cond, 0.0, u, cond_var), 1e-14)
          << "cond_var " << cond_var << ", u_cond " << u_cond << ", u " << u;

        Eigen::MatrixXd probe(1, 2);
        probe << (cond_var == 1 ? u_cond : p), (cond_var == 1 ? p : u_cond);
        EXPECT_NEAR(grid.inverse_integrate_1d(probe, cond_var)(0), u, 1e-9)
          << "cond_var " << cond_var << ", u_cond " << u_cond << ", u " << u;
      }
    }
  }
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

namespace test_tools_interpolation {
using namespace vinecopulib;
using tools_interpolation::InterpolationGrid;

namespace {
Eigen::VectorXd
uniform_knots(int m)
{
  return Eigen::VectorXd::LinSpaced(m, 0.0, 1.0);
}

Eigen::VectorXd
tail_knots(int m)
{
  // knots concentrated in the tails, as on a linear axis
  Eigen::VectorXd g(m);
  for (int i = 0; i < m; ++i) {
    const double t = static_cast<double>(i) / (m - 1);
    g(i) = 0.5 + 0.5 * std::sin((t - 0.5) * 3.14159265358979);
  }
  return g;
}
} // namespace

// A grid with different knots on its two axes interpolates a bilinear surface
// exactly, whatever the spacing of the knots.
TEST(tools_interpolation, rectangular_grid_interpolates_bilinear_exactly)
{
  auto g1 = uniform_knots(9);
  auto g2 = tail_knots(14);
  auto f = [](double x, double y) { return 1.0 + 0.5 * x - 0.3 * y + x * y; };
  Eigen::MatrixXd v(g1.size(), g2.size());
  for (int i = 0; i < g1.size(); ++i) {
    for (int j = 0; j < g2.size(); ++j) {
      v(i, j) = f(g1(i), g2(j));
    }
  }
  InterpolationGrid grid(g1, g2, v, 0);
  EXPECT_EQ(grid.get_grid_points(0), g1);
  EXPECT_EQ(grid.get_grid_points(1), g2);
  EXPECT_EQ(grid.get_values(), v);

  auto x = tools_stats::simulate_uniform(200, 2, false, { 3 });
  auto vals = grid.interpolate(x);
  for (int i = 0; i < x.rows(); ++i) {
    EXPECT_NEAR(vals(i), f(x(i, 0), x(i, 1)), 1e-12);
  }

  // flipping exchanges the axes together with their knots
  InterpolationGrid flipped = grid;
  flipped.flip();
  EXPECT_EQ(flipped.get_grid_points(0), g2);
  EXPECT_EQ(flipped.get_grid_points(1), g1);
  Eigen::MatrixXd xs = x;
  xs.col(0).swap(xs.col(1));
  EXPECT_TRUE((flipped.interpolate(xs) - vals).cwiseAbs().maxCoeff() < 1e-12);
  flipped.flip();
  EXPECT_EQ(flipped.get_values(), v);
}

// Conditional masses partition to one and agree with the conditional cdf on a
// rectangular grid, in both conditioning directions.
TEST(tools_interpolation, rectangular_grid_conditionals_are_consistent)
{
  auto g1 = uniform_knots(12);
  auto g2 = tail_knots(20);
  Eigen::MatrixXd v(g1.size(), g2.size());
  for (int i = 0; i < g1.size(); ++i) {
    for (int j = 0; j < g2.size(); ++j) {
      v(i, j) = 0.4 + std::exp(-4.0 * std::abs(g1(i) - g2(j)));
    }
  }
  InterpolationGrid grid(g1, g2, v);
  for (size_t cond_var : { 1, 2 }) {
    for (double u_cond : { 0.02, 0.37, 0.5, 0.91 }) {
      double total = 0.0;
      for (int k = 0; k < 7; ++k) {
        total +=
          grid.cond_interval_mass(u_cond, k / 7.0, (k + 1) / 7.0, cond_var);
      }
      EXPECT_NEAR(total, 1.0, 1e-13);
      Eigen::MatrixXd u(1, 2);
      const double p = 0.63;
      if (cond_var == 1) {
        u << u_cond, p;
      } else {
        u << p, u_cond;
      }
      auto q = grid.inverse_integrate_1d(u, cond_var);
      Eigen::MatrixXd uq = u;
      uq(0, cond_var == 1 ? 1 : 0) = q(0);
      EXPECT_NEAR(grid.integrate_1d(uq, cond_var)(0), p, 1e-10);
    }
  }
  // the rectangle probabilities sum to the marginal increment
  double total = 0.0;
  for (int i = 0; i < 5; ++i) {
    total += grid.rect_mass(i / 5.0, (i + 1) / 5.0, 0.2, 0.7);
  }
  EXPECT_NEAR(total, 0.5, 1e-13);
}

// On a circular axis the two end knots are the same point: when the rows at
// both ends agree, margin normalization keeps them equal, so the interpolant
// stays periodic.
TEST(tools_interpolation, normalization_keeps_periodic_ends_equal)
{
  const int m1 = 15, m2 = 11;
  auto g1 = uniform_knots(m1); // circular axis
  auto g2 = tail_knots(m2);    // linear axis
  Eigen::MatrixXd v(m1, m2);
  const double two_pi = 2 * 3.14159265358979;
  for (int i = 0; i < m1; ++i) {
    for (int j = 0; j < m2; ++j) {
      v(i, j) = 1.0 + 0.8 * std::cos(two_pi * g1(i) - 1.0) * (2 * g2(j) - 1) +
                0.3 * g2(j);
    }
  }
  ASSERT_TRUE((v.row(0) - v.row(m1 - 1)).cwiseAbs().maxCoeff() < 1e-14);
  InterpolationGrid grid(g1, g2, v);
  auto w = grid.get_values();
  EXPECT_TRUE((w.row(0) - w.row(m1 - 1)).cwiseAbs().maxCoeff() < 1e-13);
  EXPECT_TRUE((w - v).cwiseAbs().maxCoeff() > 1e-3); // it did rescale

  // the density is continuous across the cut and both margins are uniform
  Eigen::MatrixXd lo(5, 2), hi(5, 2);
  for (int k = 0; k < 5; ++k) {
    lo(k, 0) = 0.0;
    hi(k, 0) = 1.0;
    lo(k, 1) = hi(k, 1) = 0.1 + 0.2 * k;
  }
  EXPECT_TRUE(
    (grid.interpolate(lo) - grid.interpolate(hi)).cwiseAbs().maxCoeff() <
    1e-12);
  for (double u : { 0.1, 0.5, 0.9 }) {
    Eigen::MatrixXd x(1, 2);
    x << u, 1.0;
    EXPECT_NEAR(grid.integrate_2d(x)(0), u, 1e-6);
    x << 1.0, u;
    EXPECT_NEAR(grid.integrate_2d(x)(0), u, 1e-6);
  }
}
}
