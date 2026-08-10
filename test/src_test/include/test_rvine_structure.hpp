// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

namespace test_rvine_structure {

using namespace vinecopulib;

TEST(rvine_structure, triangular_array_works)
{
  TriangularArray<size_t> my_rvm({ { 2, 3, 1, 1, 1, 1 },
                                   { 1, 1, 4, 3, 2 },
                                   { 3, 2, 3, 2 },
                                   { 4, 4, 2 },
                                   { 6, 5 },
                                   { 5 } });

  EXPECT_NO_THROW(my_rvm.str());

  std::ostringstream oss;
  oss << my_rvm;
  EXPECT_EQ(oss.str(), my_rvm.str());

  std::vector<size_t> myvec = { 1, 2 };
  EXPECT_ANY_THROW(TriangularArray<size_t> my_rvm2({ { 2, 3, 1, 1, 1, 1 },
                                                     { 1, 1, 4, 3, 2 },
                                                     { 3, 2, 3, 2 },
                                                     { 4, 4, 2 },
                                                     { 6, 5 },
                                                     { 5, 2 /* TOO MANY*/ } }));
  EXPECT_ANY_THROW(TriangularArray<size_t> my_rvm3({ { 2, 3, 1, 1, 1, 1 },
                                                     { 1, 1, 4, 3, 2 },
                                                     { 3, 2, 3, 2 },
                                                     { 4, 4, 2 },
                                                     { 6, 5 },
                                                     { 5 },
                                                     { 5 } /* TOO MANY*/ }));
  EXPECT_NO_THROW(TriangularArray<size_t> my_rvm4(
    { { 2, 3, 1, 1, 1, 1 }, { 1, 1, 4, 3, 2 } }));
  my_rvm.truncate(2);
  EXPECT_EQ(my_rvm.get_trunc_lvl(), 2);

  // Different dimension
  TriangularArray<size_t> smaller_d({ { 1, 2, 3 }, { 4, 5 }, { 6 } });
  TriangularArray<size_t> larger_d(
    { { 1, 2, 3, 4 }, { 4, 5, 6 }, { 6, 7 }, { 8 } });
  EXPECT_TRUE(smaller_d < larger_d);
  EXPECT_FALSE(larger_d < smaller_d);

  // Same dimension, different truncation level
  TriangularArray<size_t> ta1({ { 1, 2, 3 }, { 4, 5 }, { 6 } });
  TriangularArray<size_t> ta2({ { 1, 2, 3 }, { 4, 5 }, { 6 } });
  ta1.truncate(2);
  ta2.truncate(3);
  EXPECT_TRUE(ta1 < ta2);
  EXPECT_FALSE(ta2 < ta1);

  // Same dimension, trunc level, and values
  TriangularArray<size_t> ta3({ { 1, 2, 3 }, { 4, 5 }, { 6 } });
  TriangularArray<size_t> ta4({ { 1, 2, 3 }, { 4, 5 }, { 6 } });
  ta3.truncate(3);
  ta4.truncate(3);
  EXPECT_FALSE(ta3 < ta4);
  EXPECT_FALSE(ta4 < ta3);
}

TEST(rvine_structure, rvine_trees_works)
{
  TriangularArray<size_t> struct_array({ { 5, 2, 6, 6, 6, 6 },
                                         { 6, 6, 1, 2, 5 },
                                         { 2, 5, 2, 5 },
                                         { 1, 1, 5 },
                                         { 3, 7 },
                                         { 7 } });
  std::vector<size_t> order = { 4, 3, 7, 1, 2, 5, 6 };

  EXPECT_NO_THROW(RVineTrees(order, struct_array));
  RVineTrees rvt(order, struct_array);

  // trees -> (order, struct_array) round-trip is the identity here
  auto rt = rvt.to_struct_array();
  EXPECT_EQ(order, rt.order);
  EXPECT_EQ(struct_array, rt.struct_array);

  // RVineStructure <-> RVineTrees round-trips
  RVineStructure rvs(order, struct_array);
  EXPECT_NO_THROW(RVineStructure{ rvt });
  EXPECT_EQ(rvs, RVineStructure(rvt));
  EXPECT_NO_THROW(RVineStructure(rvt, false));
  EXPECT_EQ(rvs, RVineStructure(rvs.get_trees())); // the get_trees() round-trip

  // truncated vine (truncate at level 3)
  TriangularArray<size_t> truncated_array(struct_array);
  truncated_array.truncate(3);
  RVineStructure rvs_trunc(order, truncated_array);

  EXPECT_NO_THROW(RVineTrees(order, truncated_array));
  RVineTrees rvt_trunc(order, truncated_array);
  auto rt_trunc = rvt_trunc.to_struct_array();
  EXPECT_EQ(order, rt_trunc.order);
  EXPECT_EQ(truncated_array, rt_trunc.struct_array);
  EXPECT_EQ(rvs_trunc, RVineStructure(rvt_trunc));
  EXPECT_EQ(rvs_trunc, RVineStructure(rvs_trunc.get_trees()));
}

TEST(rvine_structure, rvine_trees_sanity_checks)
{
  std::vector<size_t> order = { 1, 2, 3, 4 };

  // dimension mismatch between order and structure array
  TriangularArray<size_t> struct_array1({ { 2, 3 }, { 1 } });
  EXPECT_THROW(RVineTrees(order, struct_array1), std::runtime_error);

  // tree 0 does not mention every variable (variable 4 missing)
  TriangularArray<size_t> struct_array2({ { 2, 3, 3 }, { 4, 1 }, { 2 } });
  RVineTrees rvt2(order, struct_array2);
  EXPECT_THROW(rvt2.to_struct_array(), std::runtime_error);

  // proximity condition violated
  TriangularArray<size_t> struct_array3({ { 2, 3, 4 }, { 5, 1 }, { 6 } });
  RVineTrees rvt3(order, struct_array3);
  EXPECT_THROW(rvt3.to_struct_array(), std::runtime_error);
}

TEST(rvine_structure, triangular_array_conversions_work)
{
  std::vector<std::vector<size_t>> rows = {
    { 1, 2, 3, 4 }, { 4, 5, 6 }, { 6, 7 }, { 8 }
  };
  TriangularArray<size_t> ta(rows);

  // to_list() returns the rows; round-trips through the nested-vector
  // constructor.
  EXPECT_EQ(ta.to_list(), rows);
  EXPECT_TRUE(TriangularArray<size_t>(ta.to_list()) == ta);

  // to_json() carries dimension, truncation level, and rows; round-trips
  // through the JSON constructor.
  auto json = ta.to_json();
  EXPECT_EQ(json["d"].get<size_t>(), ta.get_dim());
  EXPECT_EQ(json["t"].get<size_t>(), ta.get_trunc_lvl());
  EXPECT_TRUE(TriangularArray<size_t>(json) == ta);

  // Truncated arrays keep the reduced number of rows.
  ta.truncate(2);
  EXPECT_EQ(ta.to_list().size(), static_cast<size_t>(2));
  EXPECT_EQ(ta.to_json()["t"].get<size_t>(), static_cast<size_t>(2));
  EXPECT_TRUE(TriangularArray<size_t>(ta.to_json()) == ta);

  // Arrays allocated by dimension (internal rows are over-allocated by one
  // element) still convert to rows of the logical length d - 1 - i.
  TriangularArray<size_t> filled(5, 2);
  for (size_t i = 0; i < 2; i++) {
    for (size_t j = 0; j < 5 - 1 - i; j++) {
      filled(i, j) = i + j;
    }
  }
  auto filled_rows = filled.to_list();
  ASSERT_EQ(filled_rows.size(), static_cast<size_t>(2));
  EXPECT_EQ(filled_rows[0].size(), static_cast<size_t>(4));
  EXPECT_EQ(filled_rows[1].size(), static_cast<size_t>(3));
  EXPECT_TRUE(TriangularArray<size_t>(filled_rows) == filled);
}

TEST(rvine_structure, rvine_structure_print)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;
  RVineStructure rvine_structure(mat);
  EXPECT_NO_THROW(rvine_structure.str());

  TriangularArray<size_t> my_rvm({ { 5, 2, 6, 6, 6, 6, 6 },
                                   { 6, 6, 1, 2, 5, 5 },
                                   { 2, 5, 2, 5, 2 },
                                   { 1, 1, 5, 1 },
                                   { 3, 7, 7 },
                                   { 7, 3 },
                                   { 4 } });

  std::ostringstream oss;
  oss << rvine_structure;
  EXPECT_EQ(my_rvm.str(), oss.str());
}

TEST(rvine_structure, can_convert_to_natural_order)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  TriangularArray<size_t> true_no_array({ { 6, 5, 7, 7, 7, 7 },
                                          { 7, 7, 4, 5, 6 },
                                          { 5, 6, 5, 6 },
                                          { 4, 4, 6 },
                                          { 2, 3 },
                                          { 3 } });

  RVineStructure rvine_structure(mat);
  EXPECT_EQ(rvine_structure.get_struct_array(true), true_no_array);
}

TEST(rvine_structure, min_array_is_correct)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  TriangularArray<size_t> true_min_array({ { 6, 5, 7, 7, 7, 7 },
                                           { 6, 5, 4, 5, 6 },
                                           { 5, 5, 4, 5 },
                                           { 4, 4, 4 },
                                           { 2, 3 },
                                           { 2 } });
  RVineStructure rvine_structure(mat);
  EXPECT_EQ(rvine_structure.get_min_array(), true_min_array);

  rvine_structure.truncate(2);
  EXPECT_EQ(rvine_structure.get_trunc_lvl(), 2);
  EXPECT_EQ(rvine_structure.get_min_array().get_trunc_lvl(), 2);
}

TEST(rvine_structure, needed_hfunc1_is_correct)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  TriangularArray<short unsigned> true_hfunc1({ { 0, 0, 0, 0, 1, 1 },
                                                { 0, 0, 0, 1, 1 },
                                                { 0, 0, 0, 1 },
                                                { 0, 0, 0 },
                                                { 0, 1 },
                                                { 0 } });

  RVineStructure rvine_structure(mat);
  EXPECT_EQ(rvine_structure.get_needed_hfunc1(), true_hfunc1);

  rvine_structure.truncate(2);
  EXPECT_EQ(rvine_structure.get_needed_hfunc1().get_trunc_lvl(), 2);
}

TEST(rvine_structure, needed_hfunc2_is_correct)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  TriangularArray<short unsigned> true_hfunc2({ { 1, 1, 1, 1, 1, 1 },
                                                { 1, 1, 1, 1, 1 },
                                                { 1, 1, 1, 1 },
                                                { 1, 1, 1 },
                                                { 1, 0 },
                                                { 0 } });

  RVineStructure rvine_structure(mat);
  EXPECT_EQ(rvine_structure.get_needed_hfunc2(), true_hfunc2);

  rvine_structure.truncate(2);
  EXPECT_EQ(rvine_structure.get_needed_hfunc2().get_trunc_lvl(), 2);
}

TEST(rvine_structure, construct_d_vine_struct_is_correct)
{

  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> d_vine_mat(7, 7);
  d_vine_mat << 4, 1, 5, 3, 2, 7, 7, 1, 5, 3, 2, 7, 2, 0, 5, 3, 2, 7, 3, 0, 0,
    3, 2, 7, 5, 0, 0, 0, 2, 7, 1, 0, 0, 0, 0, 7, 4, 0, 0, 0, 0, 0, 6, 0, 0, 0,
    0, 0, 0;

  std::vector<size_t> order = { 6, 4, 1, 5, 3, 2, 7 };
  RVineStructure rvine_structure(order);
  EXPECT_EQ(rvine_structure.get_matrix(), d_vine_mat);
}

TEST(rvine_structure, rvine_struct_sanity_checks_work)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  // should pass without errors
  auto rvm = RVineStructure(mat);
  auto wrong_mat = mat;

  // must be quadratic
  wrong_mat = mat.block(0, 0, 4, 5);
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

  // lower right triangle must contain zeros
  wrong_mat = mat;
  wrong_mat(6, 6) = 1;
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

  // upper left triangle must only contain 1, ..., d
  wrong_mat = mat;
  wrong_mat(0, 0) = 9;
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));
  wrong_mat(0, 0) = 0;
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

  // diagonal elements cannot appear further to the right
  wrong_mat = mat;
  wrong_mat(0, 1) = 4;
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

  // all numbers in a column most appear in each column to the left
  wrong_mat = mat;
  wrong_mat(0, 0) = 4;
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));

  // proximity condition
  wrong_mat = mat;
  wrong_mat(3, 1) = 7;
  wrong_mat(4, 1) = 1;
  EXPECT_ANY_THROW(rvm = RVineStructure(wrong_mat));
}

TEST(rvine_structure, proximity_condition_large_d)
{
  // exercise the proximity check in a higher dimension than the 7x7 case
  // above, including the error path
  const size_t d = 70;
  auto mat = CVineStructure(tools_stl::seq_int(1, d)).get_matrix();

  // a valid structure must pass the check
  EXPECT_NO_THROW(RVineStructure{ mat }); // check = true by default

  // swapping two non-antidiagonal entries within a column keeps the column's
  // value set and the diagonal intact, so only the proximity check rejects it
  // (raises the "proximity condition violated" error)
  auto wrong_mat = mat;
  std::swap(wrong_mat(0, 0), wrong_mat(2, 0));
  EXPECT_ANY_THROW(RVineStructure{ wrong_mat });
}

TEST(rvine_structure, random_sampling)
{
  for (size_t i = 0; i < 20; i++) {
    RVineStructure test =
      RVineStructure::simulate(10, false, { static_cast<int>(i + 1) });
  }
}

TEST(rvine_structure, dvine_structure)
{
  DVineStructure test(tools_stl::seq_int(1, 5));
  DVineStructure test_tr(tools_stl::seq_int(1, 5), 2);
  EXPECT_EQ(test.get_trunc_lvl(), 4);
  EXPECT_EQ(test_tr.get_trunc_lvl(), 2);
}

TEST(rvine_structure, cvine_structure)
{
  CVineStructure test(tools_stl::seq_int(1, 5));
  CVineStructure test_tr(tools_stl::seq_int(1, 5), 2);
  EXPECT_EQ(test.get_trunc_lvl(), 4);
  EXPECT_EQ(test_tr.get_trunc_lvl(), 2);
}
}
