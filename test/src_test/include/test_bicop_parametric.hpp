// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "parbicop_test.hpp"
#include "rscript.hpp"

namespace test_bicop_parametric {
using namespace vinecopulib;
using namespace tools_stl;
std::vector<int> rotations = {0, 90};

// Test if the C++ implementation of the basic methods is correct
TEST_P(ParBicopTest, parametric_bicop_is_correct) {
    if (needs_check_) {
        std::string cmd = std::string(RSCRIPT) + std::string(TEST_BICOP);
        cmd += " " + std::to_string(get_n());
        cmd += " " + std::to_string(get_family());
        cmd += " " + std::to_string(get_par());
        cmd += " " + std::to_string(get_par2());
        int sys_exit_code = system(cmd.c_str());
        if (sys_exit_code != 0) {
            throw std::runtime_error("error in system call");
        }

        int n = get_n();

        Eigen::MatrixXd results = tools_eigen::read_matxd("temp");

        // Remove temp file
        cmd = rm + "temp";
        sys_exit_code += system(cmd.c_str());
        if (sys_exit_code != 0) {
            throw std::runtime_error("error in system call");
        }

        Eigen::VectorXd par = bicop_.get_parameters();
        auto absdiff = fabs(bicop_.parameters_to_tau(par) - results(0, 0));
        ASSERT_TRUE(absdiff < 1e-4) << bicop_.str();

        // Get u-data
        Eigen::MatrixXd u = results.block(0, 1, n, 2);

        // evaluate pdf in C++
        Eigen::VectorXd f = bicop_.pdf(u);
        // assert approximate equality
        ASSERT_TRUE(f.isApprox(results.block(0, 3, n, 1), 1e-4))
                            << bicop_.str();

        // evaluate cdf in C++
        f = bicop_.cdf(u);
        // assert approximate equality
        ASSERT_TRUE(f.isApprox(results.block(0, 4, n, 1), 1e-4))
                            << bicop_.str();

        // evaluate hfunc1 in C++
        f = bicop_.hfunc1(u);
        // assert approximate equality
        ASSERT_TRUE(f.isApprox(results.block(0, 5, n, 1), 1e-4))
                            << bicop_.str();

        // evaluate hfunc2 in C++
        f = bicop_.hfunc2(u);
        // assert approximate equality
        ASSERT_TRUE(f.isApprox(results.block(0, 6, n, 1), 1e-4))
                            << bicop_.str();

        // evaluate hinv1 in C++
        f = bicop_.hinv1(u);
        // assert approximate equality
        ASSERT_TRUE(f.isApprox(results.block(0, 7, n, 1), 1e-4))
                            << bicop_.str();

        // evaluate hinv2 in C++
        f = bicop_.hinv2(u);
        // assert approximate equality
        ASSERT_TRUE(f.isApprox(results.block(0, 8, n, 1), 1e-4))
                            << bicop_.str();

        u(0, 0) = std::numeric_limits<double>::quiet_NaN();
        u(1, 1) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_NO_THROW(bicop_.pdf(u.block(0, 0, 10, 2))) << bicop_.str();
        EXPECT_TRUE(bicop_.pdf(u.block(0, 0, 1, 2)).array().isNaN()(0))
                        << bicop_.str();
        EXPECT_NO_THROW(bicop_.cdf(u.block(0, 0, 10, 2))) << bicop_.str();
        EXPECT_TRUE(bicop_.cdf(u.block(0, 0, 1, 2)).array().isNaN()(0))
                        << bicop_.str();
        EXPECT_NO_THROW(bicop_.hfunc1(u.block(0, 0, 10, 2))) << bicop_.str();
        EXPECT_TRUE(bicop_.hfunc1(u.block(0, 0, 1, 2)).array().isNaN()(0))
                        << bicop_.str();
        EXPECT_NO_THROW(bicop_.hinv1(u.block(0, 0, 10, 2))) << bicop_.str();
        EXPECT_TRUE(bicop_.hinv1(u.block(0, 0, 1, 2)).array().isNaN()(0))
                        << bicop_.str();
        EXPECT_NO_THROW(bicop_.hfunc2(u.block(0, 0, 10, 2))) << bicop_.str();
        EXPECT_TRUE(bicop_.hfunc2(u.block(0, 0, 1, 2)).array().isNaN()(0))
                        << bicop_.str();
        EXPECT_NO_THROW(bicop_.hinv2(u.block(0, 0, 10, 2))) << bicop_.str();
        EXPECT_TRUE(bicop_.hinv2(u.block(0, 0, 1, 2)).array().isNaN()(0))
                        << bicop_.str();
        EXPECT_NO_THROW(bicop_.loglik(u.block(0, 0, 10, 2))) << bicop_.str();
    }
}

// Test if the C++ implementation of select method using the mle is correct
TEST_P(ParBicopTest, bicop_select_mle_bic_is_correct) {
    std::vector<int> positive_rotations = {0, 180};
    auto true_family = bicop_.get_family_name();
    auto true_rotation = bicop_.get_rotation();
    FitControlsBicop controls({BicopFamily::indep, BicopFamily::gaussian,
                               bicop_.get_family()},
                              "mle",
                              "quadratic",
                              1.0,
                              "bic");

    if (needs_check_) {
        auto data = bicop_.simulate(get_n());
        auto bicop = Bicop(data, controls);

        //std::cout << bicop_.str() << std::endl;
        //std::cout << bicop.str() << std::endl;

        auto selected_family = bicop.get_family_name();
        EXPECT_EQ(selected_family, true_family) <<
                                                bicop_.str() << std::endl <<
                                                bicop.bic(data) << " "
                                                << bicop_.bic(data);

        if (is_member(bicop_.get_family(), bicop_families::bb)) {
            int rot_sel = bicop.get_rotation();
            if (is_member(true_rotation, positive_rotations)) {
                EXPECT_TRUE(is_member(rot_sel, positive_rotations));
            } else {
                EXPECT_FALSE(is_member(rot_sel, positive_rotations));
            }
        } else {
            EXPECT_EQ(bicop.get_rotation(), true_rotation) <<
                                                           bicop_.str()
                                                           << std::endl <<
                                                           bicop.bic(data)
                                                           << " "
                                                           << bicop_.bic(data);
        }
    }
}

// Test if the C++ implementation of select method using itau is correct
TEST_P(ParBicopTest, bicop_select_itau_bic_is_correct) {
    if (is_member(bicop_.get_family(), bicop_families::itau)) {
        auto true_family = bicop_.get_family_name();
        auto true_rotation = bicop_.get_rotation();
        FitControlsBicop controls({BicopFamily::indep, BicopFamily::gaussian,
                                   bicop_.get_family()},
                                  "itau",
                                  "quadratic",
                                  1.0,
                                  "bic");

        if (needs_check_) {
            auto data = bicop_.simulate(get_n());
            auto bicop = Bicop(data, controls);
            auto selected_family = bicop.get_family_name();
            EXPECT_EQ(selected_family, true_family) << bicop.bic(data) << " "
                                                    << bicop_.bic(data);
            EXPECT_EQ(bicop.get_rotation(), true_rotation) << bicop_.str()
                                                           << std::endl <<
                                                           bicop.bic(data)
                                                           << " "
                                                           << bicop_.bic(data);
        }
    }
}

INSTANTIATE_TEST_CASE_P(
    ParBicopTest,
    ParBicopTest,
    testing::Combine(
        testing::ValuesIn(bicop_families::parametric),
        testing::ValuesIn(rotations)
    )
);
}
