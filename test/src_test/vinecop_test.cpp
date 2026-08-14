// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/vinecop_test.hpp"

void
VinecopTest::SetUp()
{
  SKIP_WITHOUT_RSCRIPT();

  // write the exchange files for the test using VineCopula
  test_r_parity::RTempDir dir;
  dir.run(TEST_VINECOP);

  // vine structures (C++ representation reverses rows)
  model_matrix = vinecopulib::tools_eigen::read_matxs(dir.file("temp2").c_str())
                   .colwise()
                   .reverse();
  vc_matrix = vinecopulib::tools_eigen::read_matxs(dir.file("temp3").c_str())
                .colwise()
                .reverse();

  // u, pdf, sim
  Eigen::MatrixXd temp =
    vinecopulib::tools_eigen::read_matxd(dir.file("temp").c_str());
  size_t n = temp.rows();
  size_t m = model_matrix.rows();
  u = temp.block(0, 0, n, m);
  f = temp.block(0, m, n, 1);
  sim = temp.block(0, m + 1, n, m);
}
