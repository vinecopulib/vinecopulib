// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "vinecop_test.hpp"
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace test_vinecop_class {
    using namespace vinecopulib;

    TEST_F(VinecopTest, constructors_without_error) {
        Vinecop vinecop(5);
        Vinecop vinecop_indep(model_matrix);
    
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::clayton, 90);
            }
        }
    
        Vinecop vinecop_parametrized(pair_copulas, model_matrix);
    }
    
    TEST_F(VinecopTest, getters_are_correct) {
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::clayton, 90);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);
    
        for (auto& tree : vinecop.get_all_families()) {
            for (auto& fam : tree) {
                EXPECT_EQ(fam, BicopFamily::clayton);
            }
        }
    
        for (auto& tree : vinecop.get_all_pair_copulas()) {
            for (auto& pc : tree) {
                EXPECT_EQ(pc.get_family(), BicopFamily::clayton);
                EXPECT_EQ(pc.get_rotation(), 90);
            }
        }
    
        for (auto& tree : vinecop.get_all_parameters()) {
            for (auto& par : tree) {
                EXPECT_EQ(par.size(), 1);
                EXPECT_EQ(par(0), 0);
            }
        }
    
        for (auto& tree : vinecop.get_all_rotations()) {
            for (auto& rot : tree) {
                EXPECT_EQ(rot, 90);
            }
        }
    }
    
    TEST_F(VinecopTest, pdf_is_correct) {
    
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        auto par = Eigen::VectorXd::Constant(1, 3.0);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::clayton, 270, par);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);
    
        ASSERT_TRUE(vinecop.pdf(u).isApprox(f, 1e-4));
    }
    
    TEST_F(VinecopTest, cdf_is_correct) {
    
        // Create a bivariate copula and a corresponding vine with two variables
        auto pair_copulas = Vinecop::make_pair_copula_store(2);
        auto par = Eigen::VectorXd::Constant(1, 0.5);
        auto bicop = Bicop(BicopFamily::gaussian, 0, par);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::gaussian, 0, par);
            }
        }
        Eigen::Matrix<size_t, 2, 2> matrix;
        matrix << 1, 1,
                  2, 0;
        Vinecop vinecop(pair_copulas, matrix);
    
        // Test whether the analytic and simulated versions are "close" enough
        auto U = vinecop.simulate(1e1);
        ASSERT_TRUE(vinecop.cdf(U, 1e5).isApprox(bicop.cdf(U), 1e-2));
    }
    
    TEST_F(VinecopTest, simulate_is_correct) {
    
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        auto par = Eigen::VectorXd::Constant(1, 3.0);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::clayton, 270, par);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);
    
        vinecop.simulate(10);  // only check if it works
        // check the underlying transformation from independent samples
        ASSERT_TRUE(vinecop.inverse_rosenblatt(u).isApprox(sim, 1e-4));
    }
    
    TEST_F(VinecopTest, aic_bic_are_correct) {
    
        int d = 7;
        auto data = tools_stats::simulate_uniform(1000, 7);
        Vinecop true_model(d);
    
        auto pair_copulas = Vinecop::make_pair_copula_store(d);
        auto par = Eigen::VectorXd::Constant(1, 3.0);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::clayton, 0, par);
            }
        }
        Vinecop complex_model(pair_copulas, model_matrix);
    
        ASSERT_TRUE(true_model.aic(data) < complex_model.aic(data));
        ASSERT_TRUE(true_model.bic(data) < complex_model.bic(data));
    }
    
    TEST_F(VinecopTest, family_select_finds_true_rotations) {
    
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        auto par = Eigen::VectorXd::Constant(1, 3.0);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::clayton, 270, par);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);
        auto u = vinecop.simulate(10000);
        
        auto controls = FitControlsVinecop({BicopFamily::clayton}, "itau");
        // controls.set_show_trace(true);
        Vinecop fit(u, model_matrix, controls);
            
        // don't check last two trees to avoid random failures because of
        // estimation uncertainty
        auto true_rots = vinecop.get_all_rotations();
        auto fitd_rots = fit.get_all_rotations();
        true_rots.erase(true_rots.end() - 2, true_rots.end());
        fitd_rots.erase(fitd_rots.end() - 2, fitd_rots.end());
        EXPECT_EQ(true_rots, fitd_rots);
    }
    
    TEST_F(VinecopTest, select_finds_right_structure) {
        // check whether the same structure appears if we only allow for
        // independence (pair-copula estimates differ otherwise)
    
        // select structure and get matrix
        Vinecop fit(7);
        fit.select_all(u, FitControlsVinecop({BicopFamily::indep}));
        auto vcl_matrix = fit.get_matrix();
    
        // check if the same conditioned sets appear for each tree
        std::vector<std::vector<std::vector<size_t>>> vc_sets(6), vcl_sets(6);
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
                    if (tools_stl::is_same_set(s1, s2)) {
                        is_in_both = true;
                    }
                }
                if (!is_in_both) {
                    ++pairs_unequal;
                }
            }
        }
        EXPECT_EQ(pairs_unequal, 0);
    }
    
    // in what follows, we only check if funs run without error ----------
    TEST_F(VinecopTest, sparse_threshold_selection) {
        FitControlsVinecop controls(bicop_families::itau, "itau");
        controls.set_select_threshold(true);
        // controls.set_show_trace(true);
        Vinecop fit(7);
        fit.select_all(u, controls);
        fit.select_families(u, controls);
    }
    
    TEST_F(VinecopTest, sparse_truncation_selection) {
        FitControlsVinecop controls(bicop_families::itau, "itau");
        controls.set_select_truncation_level(true);
        // controls.set_show_trace(true);
        Vinecop fit(7);
        fit.select_all(u, controls);
        fit.select_families(u, controls);
    }
    
    TEST_F(VinecopTest, sparse_both_selection) {
        FitControlsVinecop controls(bicop_families::itau, "itau");
        controls.set_select_truncation_level(true);
        controls.set_select_threshold(true);
        // controls.set_show_trace(true);
        Vinecop fit(7);
        fit.select_all(u, controls);
        fit.select_families(u, controls);
    }
}
