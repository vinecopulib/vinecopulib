// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "parbicop_test.hpp"
#include "rscript.hpp"

namespace test_bicop_parametric {
    using namespace vinecopulib;
    using namespace tools_stl;

    // Test if the C++ implementation of the basic methods is correct
    TEST_P(ParBicopTest, parametric_bicop_is_correct) {
        std::string command = std::string(RSCRIPT) + "../test/test_bicop_parametric.R";
        command = command + " "  + std::to_string(get_n());
        command = command + " "  + std::to_string(get_family());
        command = command + " "  + std::to_string(get_par());
        command = command + " "  + std::to_string(get_par2());
        int sys_exit_code = system(command.c_str());
        if (sys_exit_code != 0) {
            throw std::runtime_error("error in system call");
        }

        if (needs_check_) {
            int n = get_n();

            Eigen::MatrixXd results = read_matxd("temp");
            Eigen::VectorXd par = bicop_.get_parameters();
            ASSERT_TRUE(fabs(bicop_.parameters_to_tau(par) -
                                     results(0,0)) < 1e-4);

            // evaluate pdf in C++
            Eigen::VectorXd f = bicop_.pdf(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,3,n,1), 1e-4));

            // evaluate hfunc1 in C++
            f = bicop_.hfunc1(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,4,n,1), 1e-4));

            // evaluate hfunc2 in C++
            f = bicop_.hfunc2(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,5,n,1), 1e-4));

            // evaluate hinv1 in C++
            f = bicop_.hinv1(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,6,n,1), 1e-4));

            // evaluate hinv2 in C++
            f = bicop_.hinv2(results.block(0,1,n,2));

            // assert approximate equality
            ASSERT_TRUE(f.isApprox(results.block(0,7,n,1), 1e-4));
        }
    }

    // Test if the C++ implementation of select method using the mle is correct
    TEST_P(ParBicopTest, bicop_select_mle_bic_is_correct) {
        std::vector<int> positive_rotations = {0, 180};
        std::string selection_criterion = "bic";
        auto true_family = bicop_.get_family();
        auto true_rotation = bicop_.get_rotation();
        std::vector<BicopFamily> family_set = {
                BicopFamily::indep,
                BicopFamily::gaussian,
                true_family
        };

        if (needs_check_) {
            auto data = bicop_.simulate(get_n());
            auto bicop = Bicop(data, family_set, "mle");

            auto selected_family = bicop.get_family();
            EXPECT_EQ(selected_family, true_family) <<
                            bicop.bic(data) << " " << bicop_.bic(data);

            if (is_member(true_family, bicop_families::BB)) {
                int rot_sel = bicop.get_rotation();
                if (is_member(true_rotation, positive_rotations)) {
                    EXPECT_TRUE(is_member(rot_sel, positive_rotations));
                } else {
                    EXPECT_FALSE(is_member(rot_sel, positive_rotations));
                }
            } else {
                EXPECT_EQ(bicop.get_rotation(), true_rotation) <<
                             bicop.bic(data) << " " << bicop_.bic(data);
            }
        }
    }

    // Test if the C++ implementation of select method using itau is correct
    TEST_P(ParBicopTest, bicop_select_itau_bic_is_correct) {
        if (is_member(bicop_.get_family(), bicop_families::itau)) {
            std::string selection_criterion = "bic";
            auto true_family = bicop_.get_family();
            auto true_rotation = bicop_.get_rotation();
            std::vector<BicopFamily> family_set = {
                    BicopFamily::indep,
                    BicopFamily::gaussian,
                    true_family
            };

            if (needs_check_) {
                auto data = bicop_.simulate(get_n());
                auto bicop = Bicop(data, family_set, "itau");
                auto selected_family = bicop.get_family();
                EXPECT_EQ(selected_family, true_family) << bicop.bic(data) << " " << bicop_.bic(data);
                EXPECT_EQ(bicop.get_rotation(), true_rotation) << bicop.bic(data) << " " << bicop_.bic(data);
            }
        }
    }

    INSTANTIATE_TEST_CASE_P(
            ParBicopTest,
            ParBicopTest,
            testing::Combine(
                    testing::ValuesIn(bicop_families::parametric),
                    testing::ValuesIn({0,90})
            )
    );
}
