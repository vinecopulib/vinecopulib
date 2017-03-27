// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/vinecop_test.hpp"
#include <stdexcept>

VinecopTest::VinecopTest() {
    // write temp files for the test using VineCopula
    std::string cmd = std::string(RSCRIPT) + std::string(TEST_VINECOP);
    int sys_exit_code = system(cmd.c_str());

    // vine structures (C++ representation reverses rows)
    model_matrix = tools_eigen::read_matxi("temp2").colwise().reverse();
    vc_matrix = tools_eigen::read_matxi("temp3").colwise().reverse();

    // u, pdf, sim
    Eigen::MatrixXd temp = tools_eigen::read_matxd("temp");
    int n = temp.rows();
    int m = model_matrix.rows();
    u = temp.block(0,0,n,m);
    f = temp.block(0,m,n,1);
    sim = temp.block(0,m+1,n,m);

    // remove temp files
    cmd = rm + "temp temp2 temp3";
    sys_exit_code += system(cmd.c_str());
    if (sys_exit_code != 0) {
        throw std::runtime_error("error in system call");
    }
}
