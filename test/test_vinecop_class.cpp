/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "gtest/gtest.h"
#include "include/vinecop_class.hpp"
#include "include/structselect_tools.hpp"
#include "src_test/include/vinecop_test.hpp"

namespace {
    
    VinecopTest vc_test;

    TEST(vinecop_class, constructors_without_error) {
        Vinecop vinecop_default(5);
        
        auto mat = vc_test.matrix;
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(3, VecXd::Constant(1, 3.0), 90);
            }
        }
        
        Vinecop vinecop_parametrized(pair_copulas, mat);
    }

    TEST(vinecop_class, pdf_is_correct) {
        auto mat = vc_test.matrix;
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(3, VecXd::Constant(1, 3.0), 270);
            }
        }
        Vinecop vinecop(pair_copulas, mat);
        
        auto true_pdf = Rcpp::as<VecXd>(vc_test.R.parseEval("RVinePDF(u, model)"));
        ASSERT_TRUE(vinecop.pdf(vc_test.u).isApprox(true_pdf, 1e-4));
    }

    TEST(vinecop_class, simulate_is_correct) {
        auto mat = vc_test.matrix;
        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(3, VecXd::Constant(1, 3.0), 270);
            }
        }
        Vinecop vinecop(pair_copulas, mat);
        
        auto true_sim = Rcpp::as<MatXd>(vc_test.R.parseEval("RVineSim(1000, model, U = u)"));
        ASSERT_TRUE(vinecop.simulate(1000, vc_test.u).isApprox(true_sim, 1e-4));
        
        vinecop.simulate(10);  // only check if it works
    }
            
    TEST(vinecop_class, select_finds_right_structure) {
        // check whether the same structure appears if we only allow for 
        // independence (pair-copula estimates differ otherwise)
        
        // select structure with VineCopula and vinecopulib
        Vinecop fit = Vinecop::select(vc_test.u, {0});
        vc_test.R.parseEval("fit <- RVineStructureSelect(u, familyset = 0)");
        
        // get matrices
        auto m = Rcpp::as<MatXi>(vc_test.R.parseEval("fit$Matrix"));
        auto vc_matrix = m.colwise().reverse();  // use vinecopulib orientation
        auto vcl_matrix = fit.get_matrix();
        
        // check if the same conditioned sets appear for each tree
        using namespace structselect_tools;
        std::vector<std::vector<std::vector<int>>> vc_sets(6), vcl_sets(6);
        for (int tree = 0; tree < 6; ++tree) {
            vc_sets[tree].resize(6 - tree);
            vcl_sets[tree].resize(6 - tree);
            for (int edge = 0; edge < 6 - tree; ++edge) {
                vc_sets[tree][edge].resize(2);
                vc_sets[tree][edge][0] = vc_matrix(tree, edge);
                vc_sets[tree][edge][1] = vc_matrix(5 - edge, edge);
                vcl_sets[tree][edge].resize(2);
                vcl_sets[tree][edge][0] = vc_matrix(tree, edge);
                vcl_sets[tree][edge][1] = vc_matrix(5 - edge, edge);
            }
            for (auto s1 : vc_sets[tree]) {
                bool is_in_both = false;
                for (auto s2 : vcl_sets[tree]) {
                    if (is_same_set(s1, s2))
                        is_in_both = true;
                }
                EXPECT_TRUE(is_in_both);
            }    
        }        
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
