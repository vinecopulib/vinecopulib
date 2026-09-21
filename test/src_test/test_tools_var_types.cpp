// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_var_types.hpp>

namespace test_tools_var_types {

using namespace vinecopulib;
using namespace vinecopulib::tools_var_types;

TEST(test_tools_var_types, predicates_classify_each_token)
{
  EXPECT_TRUE(is_linear(continuous()));
  EXPECT_TRUE(is_continuous(continuous()));
  EXPECT_FALSE(is_discrete(continuous()));
  EXPECT_FALSE(is_circular(continuous()));

  EXPECT_TRUE(is_discrete(discrete()));
  EXPECT_FALSE(is_continuous(discrete()));
  EXPECT_FALSE(is_linear(discrete()));
  EXPECT_FALSE(is_circular(discrete()));

  // circular is continuous but not linear
  EXPECT_TRUE(is_circular(circular()));
  EXPECT_TRUE(is_continuous(circular()));
  EXPECT_FALSE(is_linear(circular()));
  EXPECT_FALSE(is_discrete(circular()));

  EXPECT_TRUE(is_valid(continuous()));
  EXPECT_TRUE(is_valid(discrete()));
  EXPECT_TRUE(is_valid(circular()));
  EXPECT_FALSE(is_valid("x"));
  EXPECT_FALSE(is_valid(""));
}

TEST(test_tools_var_types, vector_helpers)
{
  const std::vector<std::string> cc = { "c", "c" };
  const std::vector<std::string> cd = { "c", "d" };
  const std::vector<std::string> dd = { "d", "d" };
  const std::vector<std::string> ac = { "a", "c" };
  const std::vector<std::string> ad = { "a", "d" };

  EXPECT_TRUE(all_continuous(cc));
  EXPECT_TRUE(all_continuous(ac));
  EXPECT_FALSE(all_continuous(cd));
  EXPECT_FALSE(all_continuous(ad));

  EXPECT_TRUE(all_linear(cc));
  EXPECT_FALSE(all_linear(ac));

  EXPECT_TRUE(all_discrete(dd));
  EXPECT_FALSE(all_discrete(cd));

  EXPECT_TRUE(any_circular(ac));
  EXPECT_TRUE(any_circular(ad));
  EXPECT_FALSE(any_circular(cd));

  EXPECT_EQ(count_discrete(cc), 0u);
  EXPECT_EQ(count_discrete(cd), 1u);
  EXPECT_EQ(count_discrete(dd), 2u);
  EXPECT_EQ(count_discrete(ad), 1u);

  EXPECT_EQ(all_continuous_types(3),
            (std::vector<std::string>{ "c", "c", "c" }));
}

TEST(test_tools_var_types, as_continuous_drops_discreteness_only)
{
  EXPECT_EQ(as_continuous({ "c", "d" }),
            (std::vector<std::string>{ "c", "c" }));
  EXPECT_EQ(as_continuous({ "d", "d" }),
            (std::vector<std::string>{ "c", "c" }));
  EXPECT_EQ(as_continuous({ "a", "d" }),
            (std::vector<std::string>{ "a", "c" }));
  EXPECT_EQ(as_continuous({ "a", "a" }),
            (std::vector<std::string>{ "a", "a" }));
}

TEST(test_tools_var_types, bicop_view_continuous_drops_discreteness_only)
{
  // BicopView::as_continuous() reports the discrete slot as continuous and
  // keeps the other slot unchanged.
  Bicop bc(
    BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5), { "c", "d" });
  BicopView view(bc);
  EXPECT_EQ(view.get_var_types(), (std::vector<std::string>{ "c", "d" }));
  EXPECT_EQ(view.as_continuous().get_var_types(),
            (std::vector<std::string>{ "c", "c" }));
  BicopView flipped(bc, true);
  EXPECT_EQ(flipped.get_var_types(), (std::vector<std::string>{ "d", "c" }));
  EXPECT_EQ(flipped.as_continuous().get_var_types(),
            (std::vector<std::string>{ "c", "c" }));
}

} // namespace test_tools_var_types
