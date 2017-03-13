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

#pragma once

#include "bicop.hpp"
#include "rvine_matrix.hpp"
#include <limits>

//! A class for vine copulas
class Vinecop {
public:
    Vinecop() {}
    Vinecop(int d);
    Vinecop(
        const std::vector<std::vector<BicopPtr>>& pair_copulas,
        const MatXi& matrix
    );

    static std::vector<std::vector<BicopPtr>> make_pair_copula_store(int d);
    static Vinecop select(
        const MatXd& data,
        std::vector<int> family_set = {0, 1, 2, 3, 4, 5, 6, 1001},
        std::string method = "mle",
        int truncation_level = std::numeric_limits<int>::max(),
        MatXi matrix = MatXi(0, 0),
        std::string selection_criterion = "bic",
        bool preselect_families = true,
        bool show_trace = false
    );

    BicopPtr get_pair_copula(int tree, int edge);
    int get_family(int tree, int edge);
    MatXi get_families();
    int get_rotation(int tree, int edge);
    MatXi get_rotations();
    VecXd get_parameters(int tree, int edge);
    MatXi get_matrix() {return vine_matrix_.get_matrix();}

    VecXd pdf(const MatXd& u);
    MatXd simulate(int n);
    MatXd simulate(int n, const MatXd& U);

private:
    int d_;
    RVineMatrix vine_matrix_;
    std::vector<std::vector<BicopPtr>> pair_copulas_;
};

VecXi inverse_permutation(const VecXi& order);
// reverse columns and rows of an Eigen::Matrix type object
template<typename Mat>
Mat to_upper_tri(Mat A) {return A.rowwise().reverse().colwise().reverse();}

