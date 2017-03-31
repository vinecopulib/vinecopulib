// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

// interface specfifc #defines can be set here 
// (R package does: #define INTERFACED_FROM_R)

#ifndef INTERFACED_FROM_R  // cout is not allowed for R packages
    #include <iostream>
    #define cout std::cout
#else
    #define cout Rcpp::Rcout
#endif
