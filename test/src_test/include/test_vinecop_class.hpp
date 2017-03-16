// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "vinecop_test.hpp"
#include "include/vinecop_class.hpp"
#include "include/tools_structselect.hpp"
#include "include/tools_stl.hpp"

namespace test_vinecop_class {
    using namespace vinecopulib;

    TEST_F(VinecopTest, constructors_without_error) {
        Vinecop vinecop_default(5);

        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(BicopFamily::Clayton, 90);
            }
        }

        Vinecop vinecop_parametrized(pair_copulas, model_matrix);
    }

    TEST_F(VinecopTest, pdf_is_correct) {

        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        auto par = Eigen::VectorXd::Constant(1, 3.0);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(BicopFamily::Clayton, 270, par);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);

        ASSERT_TRUE(vinecop.pdf(u).isApprox(f, 1e-4));
    }

    TEST_F(VinecopTest, simulate_is_correct) {

        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        auto par = Eigen::VectorXd::Constant(1, 3.0);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(BicopFamily::Clayton, 270, par);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);

        vinecop.simulate(10);  // only check if it works
        ASSERT_TRUE(vinecop.simulate(1000, u).isApprox(sim, 1e-4));
    }

    TEST_F(VinecopTest, select_finds_right_structure) {
// check whether the same structure appears if we only allow for
// independence (pair-copula estimates differ otherwise)

// select structure and get matrix
        Vinecop fit = Vinecop::select(u, {BicopFamily::Indep});
        auto vcl_matrix = fit.get_matrix();

// check if the same conditioned sets appear for each tree
        using namespace tools_structselect;
        std::vector<std::vector<std::vector<int>>> vc_sets(6), vcl_sets(6);
        int pairs_unequal = 0;
        for (int tree = 0; tree < 6; ++tree) {
            vc_sets[tree].resize(6 - tree);
            vcl_sets[tree].resize(6 - tree);
            for (int edge = 0; edge < 6 - tree; ++edge) {
                vc_sets[tree][edge].resize(2);
                vc_sets[tree][edge][0] = vc_matrix(tree, edge);
                vc_sets[tree][edge][1] = vc_matrix(6 - edge, edge);
                vcl_sets[tree][edge].resize(2);
                vcl_sets[tree][edge][0] = vcl_matrix(tree, edge);
                vcl_sets[tree][edge][1] = vcl_matrix(6 - edge, edge);
            }
            for (auto s1 : vc_sets[tree]) {
                bool is_in_both = false;
                for (auto s2 : vcl_sets[tree]) {
                    if (tools_stl::is_same_set(s1, s2))
                        is_in_both = true;
                }
                if (!is_in_both)
                    ++pairs_unequal;
            }
        }
        EXPECT_EQ(pairs_unequal, 0);
    }
}
