// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include "rscript.hpp"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_serialization {
    using namespace vinecopulib;

    TEST(serialization, bicop_serialization) {
        Eigen::MatrixXd mat;
        std::cout << "im here" << std::endl;
        boost::property_tree::write_json("temp", tools_serialization::matrix_to_ptree(mat));
        std::cout << "im here2" << std::endl;
        auto mat_node = tools_serialization::json_to_ptree("temp");
        std::cout << "im here3" << std::endl;
        auto mat2 = tools_serialization::ptree_to_matrix<double>(mat_node);
        std::cout << "im here4" << std::endl;

        auto pc = Bicop(BicopFamily::bb1);
        pc.to_json("temp");
        Bicop pc2("temp");

        // Remove temp file
        std::string cmd = rm + "temp";
        int sys_exit_code = system(cmd.c_str());
        if (sys_exit_code != 0) {
                throw std::runtime_error("error in system call");
        }

        EXPECT_EQ(pc.get_rotation(), pc2.get_rotation());
        EXPECT_EQ(pc.get_family_name(), pc2.get_family_name());
        ASSERT_TRUE(pc.get_parameters().isApprox(pc2.get_parameters(), 1e-4));
    }

    TEST(serialization, vinecop_serialization) {
        size_t d = 5;
        auto vc = Vinecop(d);
        vc.to_json("temp");
        auto vc2 = Vinecop("temp");

        // Remove temp file
        std::string cmd = rm + "temp";
        int sys_exit_code = system(cmd.c_str());
        if (sys_exit_code != 0) {
            throw std::runtime_error("error in system call");
        }

        for (size_t tree = 0; tree < d - 1; ++tree) {
            for (size_t edge = 0; edge < d - tree - 1; ++edge) {
                EXPECT_EQ(vc2.get_rotation(tree, edge), 0);
                EXPECT_EQ(vc2.get_family(tree, edge), BicopFamily::indep);
            }
        }
    }
    
}
