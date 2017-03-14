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
#include "src_test/include/bicop_parametric_test.hpp"

RInstance *rinstance_ptr = new RInstance;

namespace {
    // Test if the C++ implementation of select method using the mle and bic is correct
    TYPED_TEST(ParBicopTest, bicop_select_mle_bic_is_correct) {
        std::vector<int> rotations = {0, 90, 180, 270};
        std::vector<int> positive_rotations = {0,180};
        std::vector<int> bb_families = {7, 8, 9, 10};
        std::string selection_criterion = "bic";
        int true_family = this->par_bicop_.get_family();
        std::vector<int> family_set = {0,1,true_family};
        this->setup_parameters(rinstance_ptr);

        for (unsigned int j = 0; j < rotations.size(); ++j) {
            if (rotations[j] == 0 || rotations[j] == 180)
            {
                this->set_tau(rinstance_ptr, fabs(this->get_tau(rinstance_ptr)));
            } else
            {
                this->set_tau(rinstance_ptr, (-1)*fabs(this->get_tau(rinstance_ptr)));
            }
            this->set_rotation(rinstance_ptr, rotations[j]);
            this->setup_parameters(rinstance_ptr);

            if (this->needs_check_) {
                MatXd data = this->par_bicop_.simulate(this->get_n(rinstance_ptr));
                BicopPtr bicop = Bicop::select(data,
                                               selection_criterion,
                                               family_set,
                                               true,
                                               true,
                                               "mle");
                int selected_family = bicop->get_family();

                EXPECT_EQ(selected_family, true_family) << bicop->bic(data) << " " << this->par_bicop_.bic(data);
                if (is_member(true_family,bb_families))
                {
                    if (is_member(rotations[j],positive_rotations))
                    {
                        EXPECT_TRUE(is_member(bicop->get_rotation(),positive_rotations));
                    }
                    else
                    {
                        EXPECT_FALSE(is_member(bicop->get_rotation(),positive_rotations));
                    }
                } else
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
        this->setup_parameters(rinstance_ptr);

        for (unsigned int j = 0; j < rotations.size(); ++j) {
            if (rotations[j] == 0 || rotations[j] == 180)
            {
                this->set_tau(rinstance_ptr, fabs(this->get_tau(rinstance_ptr)));
            } else
            {
                this->set_tau(rinstance_ptr, (-1)*fabs(this->get_tau(rinstance_ptr)));
            }
            this->set_rotation(rinstance_ptr, rotations[j]);
            this->setup_parameters(rinstance_ptr);

            if (this->needs_check_ && is_member(this->par_bicop_.get_family(), no_itau_families) == false) {
                MatXd data = this->par_bicop_.simulate(this->get_n(rinstance_ptr));
                BicopPtr bicop = Bicop::select(data,
                                               selection_criterion,
                                               family_set,
                                               true,
                                               true,
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
