// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/vinecop_test.hpp"

VinecopTest::VinecopTest() {
    // write temp files for the test using VineCopula
    std::string command = std::string(RSCRIPT) + "../test/test_vinecop_parametric.R";
    system(command.c_str());

    // vine structures (C++ representation reverses rows)
    model_matrix = read_matxi("temp2").colwise().reverse();
    vc_matrix = read_matxi("temp3").colwise().reverse();

    // u, pdf, sim
    MatrixXd temp = read_matxd("temp");
    int n = temp.rows();
    int m = model_matrix.rows();
    u = temp.block(0,0,n,m);
    f = temp.block(0,m,n,1);
    sim = temp.block(0,m+1,n,m);

    // remove temp files
    command = rm + "temp temp2 temp3";
    system(command.c_str());
}
