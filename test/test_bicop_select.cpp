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
