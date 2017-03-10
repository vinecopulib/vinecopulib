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
        Vinecop fit = Vinecop::select(vc_test.u, {3});
        vc_test.R.parseEval("fit <- RVineStructureSelect(u, familyset = 3)");
        auto matrix = Rcpp::as<MatXi>(vc_test.R.parseEval("fit$Matrix"));
        std::cout << vc_test.matrix << std::endl << std::endl;
        std::cout << fit.get_matrix() << std::endl << std::endl;
        // TODO: compare matrices (after chaning matrix style)
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
