// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "src_test/include/test_bicop_sanity_checks.hpp"
#include "src_test/include/test_bicop_parametric.hpp"
#include "src_test/include/test_bicop_kernel.hpp"
#include "src_test/include/test_rvine_structure.hpp"
#include "src_test/include/test_serialization.hpp"
#include "src_test/include/test_tools_bobyqa.hpp"
#include "src_test/include/test_tools_stats.hpp"
#include "src_test/include/test_vinecop_class.hpp"
#include "src_test/include/test_vinecop_sanity_checks.hpp"
#include "src_test/include/test_weights.hpp"

using namespace test_bicop_sanity_checks;
using namespace test_bicop_parametric;
using namespace test_bicop_kernel;
using namespace test_rvine_structure;
using namespace test_serialization;
using namespace test_tools_bobyqa;
using namespace test_tools_stats;
using namespace test_vinecop_class;
using namespace test_vinecop_sanity_checks;
using namespace test_weights;

int main(int argc, char **argv) {
    // ::testing::InitGoogleTest(&argc, argv);
    // return RUN_ALL_TESTS();

    auto bc = Bicop(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 3));
    auto u = bc.simulate(1000);

    Eigen::MatrixXd unew(1000, 4);
    unew.block(0, 0, 1000, 2) = u;
    unew.block(0, 2, 1000, 2) = u;
    unew.col(2) = (unew.col(2).array() * 2).floor() / 2;
    // unew.col(3) = (unew.col(3).array() * 2).floor() / 2.0;
    bc.set_discrete_vars({0});
    bc.fit(unew);
    std::cout << bc.str() << std::endl;

    return 0;
}
