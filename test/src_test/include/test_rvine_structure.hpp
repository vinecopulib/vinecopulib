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
  EXPECT_NO_THROW(rvt.to_struct_array());

  auto my_rvt = rvt.to_struct_array();
  EXPECT_EQ(order, std::get<0>(my_rvt));
  EXPECT_EQ(struct_array, std::get<1>(my_rvt));

  RVineStructure rvs(order, struct_array);
  EXPECT_NO_THROW(RVineStructure rvs2(rvt));
  EXPECT_EQ(rvs, RVineStructure(rvt));
  EXPECT_NO_THROW(RVineStructure rvs3(rvt, false));

  // Truncated vine test (truncate at level 3)
  size_t trunc_lvl = 3;
  TriangularArray<size_t> truncated_array(struct_array);
  truncated_array.truncate(trunc_lvl);
  RVineStructure rvs_trunc(order, truncated_array);

  EXPECT_NO_THROW(RVineTrees(order, truncated_array));
  RVineTrees rvt_trunc(order, truncated_array);
  auto my_rvt_trunc = rvt_trunc.to_struct_array();
  EXPECT_EQ(order, std::get<0>(my_rvt_trunc));
  EXPECT_EQ(truncated_array, std::get<1>(my_rvt_trunc));
  EXPECT_NO_THROW(RVineStructure rvs_trunc2(rvt_trunc));
  EXPECT_EQ(rvs_trunc, RVineStructure(rvt_trunc));
}

TEST(rvine_structure, rvine_trees_sanity_checks)
{
  // throws_on_dimension_mismatch)
  std::vector<size_t> order = { 1, 2, 3, 4 };
  TriangularArray<size_t> struct_array1({ { 2, 3 }, { 1 } });

  EXPECT_THROW(RVineTrees(order, struct_array1), std::runtime_error);

  // throws_on_missing_vars_in_first_tree
  TriangularArray<size_t> struct_array2(
    { { 2, 3, 3 }, // edge (1,2), (2,3), (3,3) -- missing variable 4 if d=4
      { 4, 1 },
      { 2 } });

  RVineTrees rvt2(order, struct_array2);
  EXPECT_THROW(rvt2.to_struct_array(), std::runtime_error);

  // throws_on_proximity_condition)
  TriangularArray<size_t> struct_array3({ { 2, 3, 4 }, { 5, 1 }, { 6 } });

  RVineTrees rvt3(order, struct_array3);
  EXPECT_THROW(rvt3.to_struct_array(), std::runtime_error);
}

TEST(rvine_structure, rvine_trees_merge_constructor_works)
{
  // Define two simple vines on variables {1, 2, 3, 4} and {5, 6, 7, 8}
  TriangularArray<size_t> arr1({ { 2, 3, 4 }, { 3, 4 }, { 4 } });
  std::vector<size_t> order1 = { 1, 2, 3, 4 };

  TriangularArray<size_t> arr2({ { 6, 7, 8 }, { 7, 8 }, { 8 } });
  std::vector<size_t> order2 = { 5, 6, 7, 8 };

  RVineTrees vine1(order1, arr1);
  RVineTrees vine2(order2, arr2);

  std::vector<RVineTrees> vines = { vine1, vine2 };

  // Should not throw
  EXPECT_NO_THROW(RVineTrees merged_tmp(vines));

  RVineTrees merged(vines);

  // The merged structure has twice as many edges per tree (2x each tree)
  auto merged_trees = merged.get_trees();
  for (size_t t = 0; t < merged_trees.size(); ++t) {
    EXPECT_EQ(merged_trees[t].size(), 2 * vine1.get_trees()[t].size());
  }
  EXPECT_THROW(merged.to_struct_array(),
               std::runtime_error); // not a valid structure

  EXPECT_NO_THROW(
    merged.to_struct_array(true)); // valid structure with fill_missing = true

  auto merged_struct = merged.to_struct_array(true);
  auto merged_order = std::get<0>(merged_struct);
  auto merged_struct_array = std::get<1>(merged_struct);
  EXPECT_NO_THROW(RVineStructure(merged_order, merged_struct_array));

  std::vector<RVineTrees> empty_list;
  EXPECT_THROW(RVineTrees empty_tmp(empty_list), std::runtime_error);

  RVineTrees empty_vine;
  vines.push_back(empty_vine);
  EXPECT_THROW(RVineTrees empty_tmp2(vines), std::runtime_error);
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

TEST(rvine_structure, random_sampling)
{
  for (size_t i = 0; i < 20; i++) {
    RVineStructure test = RVineStructure::simulate(10);
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
