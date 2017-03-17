// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "parbicop_test.hpp"

namespace test_bicop_select {
    using namespace vinecopulib;
    using namespace tools_stl;
    
    // Test if the C++ implementation of select method using the mle and bic is correct
    TYPED_TEST(ParBicopTest, bicop_select_mle_bic_is_correct) {
        std::vector<int> rotations = {0, 90, 180, 270};
        std::vector<int> positive_rotations = {0, 180};
        std::string selection_criterion = "bic";
        auto true_family = this->par_bicop_.get_family();
        std::vector<BicopFamily> family_set = {
            BicopFamily::indep, 
            BicopFamily::gaussian, 
            true_family
        };

        for (auto rotation : rotations) {
            this->setup_parameters(rotation);

            if (this->needs_check_) {
                auto data = this->par_bicop_.simulate(this->get_n());
                auto bicop = Bicop::select(data, family_set, "mle");

                auto selected_family = bicop->get_family();
                EXPECT_EQ(selected_family, true_family) << 
                    bicop->bic(data) << " " << this->par_bicop_.bic(data);

                if (is_member(true_family, bicop_families::BB)) {
                    int rot_sel = bicop->get_rotation();
                    if (is_member(rotation, positive_rotations)) {
                        EXPECT_TRUE(is_member(rot_sel, positive_rotations));
                    } else {
                        EXPECT_FALSE(is_member(rot_sel, positive_rotations));
                    }
                } else {
                    EXPECT_EQ(bicop->get_rotation(), rotation) << 
                        bicop->bic(data) << " " << this->par_bicop_.bic(data);
                }
            }
        }
    }

    TYPED_TEST(ParBicopTest, bicop_select_itau_bic_is_correct) {
        if (is_member(this->par_bicop_.get_family(), bicop_families::itau)) {
            std::vector<int> no_itau_families = {7, 8, 9, 10};
            std::vector<int> rotations = {0, 90, 180, 270};
            std::string selection_criterion = "bic";
            auto true_family = this->par_bicop_.get_family();
            std::vector<BicopFamily> family_set = {
                BicopFamily::indep, 
                BicopFamily::gaussian, 
                true_family
            };
            this->setup_parameters();

            for (auto rotation : rotations) {
                this->setup_parameters(rotation);
                if (this->needs_check_) {
                    auto data = this->par_bicop_.simulate(this->get_n());
                    auto bicop = Bicop::select(data, family_set, "itau");

                    auto selected_family = bicop->get_family();
                    EXPECT_EQ(selected_family, true_family) << 
                        bicop->bic(data) << " " << this->par_bicop_.bic(data);
                    EXPECT_EQ(bicop->get_rotation(), rotation) << 
                        bicop->bic(data) << " " << this->par_bicop_.bic(data);
                }
            }
        }
    }
}
