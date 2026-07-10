// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "test_utils.hpp"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include <vinecopulib.hpp>

namespace test_serialization {
using namespace vinecopulib;
using test_utils::all_close;

class TemporaryFile
{
public:
  explicit TemporaryFile(const std::string& filename)
    : filename_(filename)
  {
  }

  ~TemporaryFile() { std::remove(filename_.c_str()); }

  const std::string& filename() const { return filename_; }

private:
  std::string filename_;
};

inline void
expect_bicops_equal(const Bicop& expected, const Bicop& actual)
{
  EXPECT_EQ(expected.get_family(), actual.get_family());
  EXPECT_EQ(expected.get_rotation(), actual.get_rotation());
  EXPECT_EQ(expected.get_var_types(), actual.get_var_types());
  EXPECT_EQ(expected.get_npars(), actual.get_npars());
  EXPECT_TRUE(all_close(expected.get_parameters(), actual.get_parameters()));
}

inline void
expect_vinecops_equal(const Vinecop& expected, const Vinecop& actual)
{
  EXPECT_EQ(expected.get_dim(), actual.get_dim());
  EXPECT_EQ(expected.get_trunc_lvl(), actual.get_trunc_lvl());
  EXPECT_EQ(expected.get_matrix(), actual.get_matrix());
  EXPECT_EQ(expected.get_all_families(), actual.get_all_families());
  EXPECT_EQ(expected.get_all_rotations(), actual.get_all_rotations());
  EXPECT_EQ(expected.get_var_types(), actual.get_var_types());

  const auto expected_parameters = expected.get_all_parameters();
  const auto actual_parameters = actual.get_all_parameters();
  ASSERT_EQ(expected_parameters.size(), actual_parameters.size());
  for (size_t tree = 0; tree < expected_parameters.size(); ++tree) {
    ASSERT_EQ(expected_parameters[tree].size(), actual_parameters[tree].size());
    for (size_t edge = 0; edge < expected_parameters[tree].size(); ++edge) {
      EXPECT_TRUE(all_close(expected_parameters[tree][edge],
                            actual_parameters[tree][edge]));
    }
  }
}

inline Vinecop
make_mixed_vinecop()
{
  auto pair_copulas = Vinecop::make_pair_copula_store(4);

  pair_copulas[0][0] =
    Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.4));
  pair_copulas[0][1] =
    Bicop(BicopFamily::clayton, 90, Eigen::VectorXd::Constant(1, 1.2));
  pair_copulas[0][2] = Bicop(BicopFamily::tll);

  Eigen::VectorXd student_parameters(2);
  student_parameters << -0.3, 5.0;
  pair_copulas[1][0] = Bicop(BicopFamily::student, 0, student_parameters);
  pair_copulas[1][1] =
    Bicop(BicopFamily::frank, 0, Eigen::VectorXd::Constant(1, 2.0));

  Eigen::VectorXd bb1_parameters(2);
  bb1_parameters << 2.0, 1.5;
  pair_copulas[2][0] = Bicop(BicopFamily::bb1, 0, bb1_parameters);

  return Vinecop(DVineStructure({ 1, 2, 3, 4 }), pair_copulas);
}

TEST(serialization, parametric_bicop_cbor_roundtrip)
{
  Eigen::VectorXd parameters(2);
  parameters << 0.5, 4.0;
  const Bicop bicop(
    BicopFamily::student, 0, parameters, std::vector<std::string>{ "c", "d" });

  TemporaryFile automatic_file("serialization_parametric.cbor");
  bicop.to_file(automatic_file.filename());
  expect_bicops_equal(bicop, Bicop(automatic_file.filename()));

  TemporaryFile explicit_file("serialization_parametric.bin");
  bicop.to_file(explicit_file.filename(), SerializationFormat::cbor);
  expect_bicops_equal(
    bicop, Bicop(explicit_file.filename(), SerializationFormat::cbor));
}

TEST(serialization, tll_bicop_cbor_roundtrip)
{
  const auto data = tools_stats::simulate_uniform(30, 2, true, { 1 });
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_grid_size(10);
  const Bicop bicop(data, controls);

  TemporaryFile file("serialization_tll.cbor");
  bicop.to_file(file.filename());
  expect_bicops_equal(bicop, Bicop(file.filename()));
}

TEST(serialization, mixed_vinecop_cbor_roundtrip)
{
  const auto vinecop = make_mixed_vinecop();

  TemporaryFile automatic_file("serialization_vinecop.cbor");
  vinecop.to_file(automatic_file.filename());
  expect_vinecops_equal(vinecop, Vinecop(automatic_file.filename()));

  auto truncated = vinecop;
  truncated.truncate(2);
  TemporaryFile explicit_file("serialization_vinecop_truncated.bin");
  truncated.to_file(explicit_file.filename(), SerializationFormat::cbor);
  expect_vinecops_equal(
    truncated, Vinecop(explicit_file.filename(), SerializationFormat::cbor));
}

TEST(serialization, rvine_structure_cbor_roundtrip)
{
  const RVineStructure structure(DVineStructure({ 1, 2, 3, 4 }, 2));
  TemporaryFile file("serialization_structure.cbor");
  structure.to_file(file.filename());
  const RVineStructure restored(file.filename());

  EXPECT_EQ(structure.get_trunc_lvl(), restored.get_trunc_lvl());
  EXPECT_EQ(structure.get_matrix(), restored.get_matrix());
}

TEST(serialization, json_behavior_is_unchanged)
{
  const Bicop bicop(
    BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));
  TemporaryFile file("serialization_compatibility.json");
  bicop.to_file(file.filename());

  std::ifstream input(file.filename());
  const std::string contents((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  EXPECT_EQ(bicop.to_json().dump() + "\n", contents);
  expect_bicops_equal(bicop, Bicop(file.filename()));
}

TEST(serialization, malformed_files_report_the_encoding)
{
  TemporaryFile cbor_file("serialization_invalid.cbor");
  {
    std::ofstream output(cbor_file.filename(), std::ios::binary);
    output.put(static_cast<char>(0xff));
  }
  try {
    Bicop bicop(cbor_file.filename());
    FAIL() << "invalid CBOR was accepted";
  } catch (const std::runtime_error& exception) {
    EXPECT_NE(std::string(exception.what()).find("CBOR"), std::string::npos);
  }

  TemporaryFile json_file("serialization_invalid.json");
  {
    std::ofstream output(json_file.filename());
    output << "{";
  }
  try {
    Bicop bicop(json_file.filename());
    FAIL() << "invalid JSON was accepted";
  } catch (const std::runtime_error& exception) {
    EXPECT_NE(std::string(exception.what()).find("JSON"), std::string::npos);
  }
}

TEST(serialization, reports_open_and_format_errors)
{
  EXPECT_THROW(Bicop(std::string("serialization_missing/model.json")),
               std::runtime_error);

  const Bicop bicop;
  EXPECT_THROW(bicop.to_file("serialization_missing/model.cbor"),
               std::runtime_error);

  TemporaryFile file("serialization_unsupported");
  const auto unsupported = static_cast<SerializationFormat>(100);
  EXPECT_THROW(bicop.to_file(file.filename(), unsupported),
               std::invalid_argument);
}
}
