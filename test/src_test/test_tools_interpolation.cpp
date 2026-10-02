// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <string>
#include <utility>
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

//! An equally spaced grid on [0, 1] and its trapezoid weights.
std::pair<Eigen::VectorXd, Eigen::VectorXd>
grid_and_weights(int m)
{
  Eigen::VectorXd g(m);
  for (int i = 0; i < m; ++i) {
    g(i) = static_cast<double>(i) / (m - 1);
  }
  Eigen::VectorXd w(m);
  w(0) = (g(1) - g(0)) / 2;
  w(m - 1) = (g(m - 1) - g(m - 2)) / 2;
  for (int i = 1; i < m - 1; ++i) {
    w(i) = (g(i + 1) - g(i - 1)) / 2;
  }
  return { g, w };
}

//! `exp(-concentration * |g_i - g_j|)`, tilted off uniform margins.
Eigen::MatrixXd
concentrated_surface(const Eigen::VectorXd& g, double concentration)
{
  const auto m = g.size();
  Eigen::MatrixXd v(m, m);
  for (ptrdiff_t i = 0; i < m; ++i) {
    for (ptrdiff_t j = 0; j < m; ++j) {
      v(i, j) = std::exp(-concentration * std::abs(g(i) - g(j))) + 1e-3 * g(i);
    }
  }
  return v;
}

//! Both margins of the normalized `v` are uniform, and `v` and its transpose
//! normalize to transposes of each other, bit for bit.
void
expect_normalized_and_equivariant(const Eigen::VectorXd& g,
                                  const Eigen::VectorXd& w,
                                  const Eigen::MatrixXd& v,
                                  const std::string& label)
{
  const Eigen::MatrixXd normalized = InterpolationGrid(g, v).get_values();
  const double rows = ((normalized * w).array() - 1.0).abs().maxCoeff();
  const double cols =
    ((normalized.transpose() * w).array() - 1.0).abs().maxCoeff();
  EXPECT_LT(std::max(rows, cols), 1e-13) << label;
  const Eigen::MatrixXd flipped =
    InterpolationGrid(g, Eigen::MatrixXd(v.transpose())).get_values();
  EXPECT_TRUE((flipped.transpose().array() == normalized.array()).all())
    << label;
}

} // namespace

// margins converge on concentrated surfaces, and transposes normalize alike
TEST(tools_interpolation, normalization_converges_on_a_concentrated_surface)
{
  const auto [g, w] = grid_and_weights(30);
  for (double concentration : { 40.0, 200.0, 400.0 }) {
    expect_normalized_and_equivariant(g,
                                      w,
                                      concentrated_surface(g, concentration),
                                      std::to_string(concentration));
  }
}

// a support in two blocks, or with a corner of its own, converges too
TEST(tools_interpolation, normalization_converges_on_a_disconnected_surface)
{
  const auto [g, w] = grid_and_weights(30);
  for (double concentration : { 40.0, 200.0, 400.0 }) {
    Eigen::MatrixXd halves = concentrated_surface(g, concentration);
    halves.topRightCorner(15, 15).setZero();
    halves.bottomLeftCorner(15, 15).setZero();
    expect_normalized_and_equivariant(
      g, w, halves, "halves, " + std::to_string(concentration));
    Eigen::MatrixXd corner = concentrated_surface(g, concentration);
    corner.row(0).tail(29).setZero();
    corner.col(0).tail(29).setZero();
    expect_normalized_and_equivariant(
      g, w, corner, "corner, " + std::to_string(concentration));
  }
}

// the masses over a partition of the first argument add up to the strip's
TEST(tools_interpolation, rect_mass_telescopes_in_the_first_argument)
{
  auto grid = skewed_grid(30);
  for (int k : { 1, 2, 5, 16, 128, 1000 }) {
    for (double b2 : { 0.125, 0.5, 0.875, 1.0 }) {
      for (double width : { 1.0 / 8, 1.0 / 512 }) {
        const double a2 = std::max(b2 - width, 0.0);
        const double whole = grid.rect_mass(0.0, 1.0, a2, b2);
        double total = 0.0;
        for (int i = 0; i < k; ++i) {
          total += grid.rect_mass(
            static_cast<double>(i) / k, static_cast<double>(i + 1) / k, a2, b2);
        }
        EXPECT_NEAR(total, whole, 1e-14 * whole)
          << "k = " << k << ", (a2, b2) = (" << a2 << ", " << b2 << ")";
      }
    }
  }
}

// a grid and its transpose give the same masses and cdf, bit for bit
TEST(tools_interpolation, mass_and_cdf_commute_with_transposition)
{
  const int m = 30;
  auto grid = skewed_grid(m);
  Eigen::VectorXd g(m);
  for (int i = 0; i < m; ++i) {
    g(i) = static_cast<double>(i) / (m - 1);
  }
  auto flipped =
    InterpolationGrid(g, Eigen::MatrixXd(grid.get_values().transpose()), 0);
  for (double a1 : { 0.0, 0.1, 0.4 }) {
    for (double a2 : { 0.0, 0.15, 0.55 }) {
      for (double width : { 0.3, 1.0 / 512 }) {
        const double b1 = std::min(a1 + width, 1.0);
        const double b2 = std::min(a2 + 2.0 * width, 1.0);
        EXPECT_EQ(grid.rect_mass(a1, b1, a2, b2),
                  flipped.rect_mass(a2, b2, a1, b1))
          << "(a1, a2) = (" << a1 << ", " << a2 << "), width " << width;
      }
    }
  }
  Eigen::MatrixXd u(4, 2);
  u << 0.3, 0.7, 0.05, 0.95, 1.0, 0.4, 0.5, 0.5;
  const Eigen::MatrixXd swapped = u.rowwise().reverse();
  EXPECT_TRUE(
    (grid.integrate_2d(u).array() == flipped.integrate_2d(swapped).array())
      .all());
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
