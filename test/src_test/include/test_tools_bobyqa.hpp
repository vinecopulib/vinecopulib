// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_bobyqa.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <cmath>

namespace test_tools_bobyqa {

using namespace vinecopulib;

TEST(test_tools_bobyqa, const_function) {

    auto f = [](long /*n*/, const double* /*x*/) -> double {
        return 0.0;
    };

    const long variables_count = 2;
    const long number_of_interpolation_conditions = variables_count + 2;
    const long ws_size = (number_of_interpolation_conditions + 5) *
                         (number_of_interpolation_conditions
                          + variables_count) +
                         3 * variables_count * (variables_count + 5) / 2;


    double *lb = new double[variables_count];
    double *ub = new double[variables_count];
    double *x = new double[variables_count];
    double *working_space = new double[ws_size];

    lb[0] = -1.0;
    lb[1] = -1.0;
    ub[0] = 1.0;
    ub[1] = 1.0;
    x[0] = 0.0;
    x[1] = 0.0;

    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 5;

    double result = tools_bobyqa::impl(f, variables_count,
                                       number_of_interpolation_conditions,
                                       x, lb, ub,
                                       initial_trust_region_radius,
                                       final_trust_region_radius,
                                       max_function_calls_count,
                                       working_space);

    ASSERT_TRUE(fabs(result) < 1e-6);
    ASSERT_TRUE(fabs(x[0]) < 1e-6);
    ASSERT_TRUE(fabs(x[1]) < 1e-6);

    // delete dynamically allocated objects
    delete[] x;
    delete[] lb;
    delete[] ub;
    delete[] working_space;
}

TEST(test_tools_bobyqa, complex_quadratic_function) {

    auto f = [](long /*n*/, const double *x) -> double {
        return -4*x[0]*x[1] + 5*x[0]*x[0] + 8*x[1]*x[1] +
            16*sqrt(5.0)*x[0] + 8*sqrt(5.0)*x[1] - 44.0;
    };

    const long variables_count = 2;
    const long number_of_interpolation_conditions =
        (variables_count + 1)*(variables_count + 2)/2;
    const long ws_size = (number_of_interpolation_conditions + 5) *
                          (number_of_interpolation_conditions
                           + variables_count) +
                          3 * variables_count * (variables_count + 5) / 2;


    double *lb = new double[variables_count];
    double *ub = new double[variables_count];
    double *x = new double[variables_count];
    double *working_space = new double[ws_size];

    lb[0] = -4.0;
    lb[1] = -3.0;
    ub[0] = 5.0;
    ub[1] = 5.0;
    x[0] = 0.0;
    x[1] = -sqrt(5.0);

    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 22;

    double result = tools_bobyqa::impl(f, variables_count,
                                       number_of_interpolation_conditions,
                                       x, lb, ub,
                                       initial_trust_region_radius,
                                       final_trust_region_radius,
                                       max_function_calls_count,
                                       working_space);

    ASSERT_TRUE(fabs(result + 142.9968943799848) < 1e-6);
    ASSERT_TRUE(fabs(x[0] + 4.0) < 1e-6);
    ASSERT_TRUE(fabs(x[1] + 2.11803398875050552) < 1e-6);

    // delete dynamically allocated objects
    delete[] x;
    delete[] lb;
    delete[] ub;
    delete[] working_space;
}


TEST(test_tools_bobyqa, quadratic_function_with_jump) {

    auto f = [](long /*n*/, const double *x) -> double {
        return x[1] < 0.5 ? x[0] * x[0] +
            x[1] * x[1] : x[0] * x[0] + x[1] * x[1] + 10;
    };

    const long variables_count = 2;
    const long number_of_interpolation_conditions =
        (variables_count + 1)*(variables_count + 2)/2;
    const long ws_size = (number_of_interpolation_conditions + 5) *
                         (number_of_interpolation_conditions
                          + variables_count) +
                         3 * variables_count * (variables_count + 5) / 2;


    double *lb = new double[variables_count];
    double *ub = new double[variables_count];
    double *x = new double[variables_count];
    double *working_space = new double[ws_size];

    lb[0] = -1.0;
    lb[1] = -1.0;
    ub[0] = 1.0;
    ub[1] = 1.0;
    x[0] = 0.5;
    x[1] = 0.5;

    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 8;

    double result = tools_bobyqa::impl(f, variables_count,
                                       number_of_interpolation_conditions,
                                       x, lb, ub,
                                       initial_trust_region_radius,
                                       final_trust_region_radius,
                                       max_function_calls_count,
                                       working_space);

    ASSERT_TRUE(fabs(result - 0.497004933606359666) < 1e-6);
    ASSERT_TRUE(fabs(x[0] - 0.49899993347109322661) < 1e-6);
    ASSERT_TRUE(fabs(x[1] - 0.49800000000221306128) < 1e-6);

    // delete dynamically allocated objects
    delete[] x;
    delete[] lb;
    delete[] ub;
    delete[] working_space;
}

}
