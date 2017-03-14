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

#include "include/vinecop_test.hpp"

VinecopTest::VinecopTest() {
    // write temp files for the test using VineCopula
    std::string command = std::string(RSCRIPT) + "../test/test_vinecop_parametric.R";
    system(command.c_str());

    // vine structures (C++ representation reverses rows)
    model_matrix = read_matxi("temp2").colwise().reverse();
    vc_matrix = read_matxi("temp3").colwise().reverse();

    // u, pdf, sim
    MatXd temp = read_matxd("temp");
    int n = temp.rows();
    int m = model_matrix.rows();
    u = temp.block(0,0,n,m);
    f = temp.block(0,m,n,1);
    sim = temp.block(0,m+1,n,m);

    // remove temp files
    system("rm temp temp2 temp3");
}