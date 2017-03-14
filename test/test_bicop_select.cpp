// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "src_test/include/parbicop_test.hpp"

namespace {
    // Test if the C++ implementation of select method using the mle and bic is correct
    TYPED_TEST(ParBicopTest, bicop_select_mle_bic_is_correct) {
        std::vector<int> rotations = {0, 90, 180, 270};
        std::vector<int> positive_rotations = {0,180};
        std::vector<int> bb_families = {7, 8, 9, 10};
        std::string selection_criterion = "bic";
        int true_family = this->par_bicop_.get_family();
        std::vector<int> family_set = {0,1,true_family};

        for (unsigned int j = 0; j < rotations.size(); ++j) {
            this->setup_parameters(rotations[j]);

            if (this->needs_check_)
            {
                MatXd data = this->par_bicop_.simulate(this->get_n());
                BicopPtr bicop = Bicop::select(data,
                                               family_set,
                                               "mle");

                int selected_family = bicop->get_family();
                EXPECT_EQ(selected_family, true_family) << bicop->bic(data) << " " << this->par_bicop_.bic(data);

                if (tools_stl::is_member(true_family,bb_families))
                {
                    if (tools_stl::is_member(rotations[j],positive_rotations))
                    {
                        EXPECT_TRUE(tools_stl::is_member(bicop->get_rotation(),positive_rotations));
                    }
                    else
                    {
                        EXPECT_FALSE(tools_stl::is_member(bicop->get_rotation(),positive_rotations));
                    }
                }
                else
                {
                    EXPECT_EQ(bicop->get_rotation(), rotations[j]) << bicop->bic(data) << " " << this->par_bicop_.bic(data);
                }
            }
        }
    }

    TYPED_TEST(ParBicopTest, bicop_select_itau_bic_is_correct) {
        std::vector<int> no_itau_families = {7, 8, 9, 10};
        std::vector<int> rotations = {0, 90, 180, 270};
        std::string selection_criterion = "bic";
        int true_family = this->par_bicop_.get_family();
        std::vector<int> family_set = {0,1,true_family};
        this->setup_parameters();

        for (unsigned int j = 0; j < rotations.size(); ++j) {
            this->setup_parameters(rotations[j]);

            if (this->needs_check_ && tools_stl::is_member(this->par_bicop_.get_family(), no_itau_families) == false) {
                MatXd data = this->par_bicop_.simulate(this->get_n());
                BicopPtr bicop = Bicop::select(data,
                                               family_set,
                                               "itau");

                int selected_family = bicop->get_family();
                EXPECT_EQ(selected_family, true_family) << bicop->bic(data) << " " << this->par_bicop_.bic(data);
                EXPECT_EQ(bicop->get_rotation(), rotations[j]) << bicop->bic(data) << " " << this->par_bicop_.bic(data);
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
