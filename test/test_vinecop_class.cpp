/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
*
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "include/vinecop_class.hpp"
#include "include/structselect_tools.hpp"
#include "src_test/include/vinecop_test.hpp"

namespace {

    TEST_F(VinecopTest, constructors_without_error) {
        Vinecop vinecop_default(5);

        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(3, VecXd::Constant(1, 3.0), 90);
            }
        }

        Vinecop vinecop_parametrized(pair_copulas, model_matrix);
    }

    TEST_F(VinecopTest, pdf_is_correct) {

        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(3, VecXd::Constant(1, 3.0), 270);
            }
        }
        Vinecop vinecop(pair_copulas, model_matrix);

        ASSERT_TRUE(vinecop.pdf(u).isApprox(f, 1e-4));
    }

    TEST_F(VinecopTest, simulate_is_correct) {

        auto pair_copulas = Vinecop::make_pair_copula_store(7);
        for (auto& tree : pair_copulas) {
            for (auto& pc : tree) {
                pc = Bicop::create(3, VecXd::Constant(1, 3.0), 270);
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
        Vinecop fit = Vinecop::select(u, {0});
        auto vcl_matrix = fit.get_matrix();

        std::cout << vc_matrix << std::endl << std::endl;
        std::cout << vcl_matrix << std::endl << std::endl;

        // check if the same conditioned sets appear for each tree
        using namespace structselect_tools;
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
                    if (is_same_set(s1, s2))
                        is_in_both = true;
                }
                if (!is_in_both)
                    ++pairs_unequal;
            }
        }
        EXPECT_EQ(pairs_unequal, 0);
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
