// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/distribution_tools.hpp"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>

template <typename T>
void time(const std::string label, const T &it)
{
    auto start = std::chrono::high_resolution_clock::now();
    it();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << label << "Elapsed time: " << elapsed.count() << " s\n";
}


int main(int argc, char **argv) {
    Eigen::MatrixXd m = Eigen::MatrixXd::Ones(10000, 10000);
    m *= 0.5;

    time("Boost dnorm:", [m]{ auto a = dnorm(m); });
    time("Boost pnom:", [m]{ auto a = pnorm(m); });
    time("Boost qnom:", [m]{ auto a = qnorm(m); });

    return 0;
}
