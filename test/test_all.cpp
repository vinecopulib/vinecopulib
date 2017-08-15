// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "src_test/include/test_bicop_sanity_checks.hpp"
#include "src_test/include/test_bicop_parametric.hpp"
#include "src_test/include/test_bicop_kernel.hpp"
#include "src_test/include/test_rvine_matrix.hpp"
#include "src_test/include/test_serialization.hpp"
#include "src_test/include/test_tools_stats.hpp"
#include "src_test/include/test_vinecop_class.hpp"
#include "src_test/include/test_vinecop_sanity_checks.hpp"

using namespace test_bicop_sanity_checks;
using namespace test_bicop_parametric;
using namespace test_bicop_kernel;
using namespace test_rvine_matrix;
using namespace test_serialization;
using namespace test_tools_stats;
using namespace test_vinecop_class;
using namespace test_vinecop_sanity_checks;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
