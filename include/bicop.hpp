// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

// differentiate between compilation of the library and of the client code
#if defined (_WIN32) 
#    if defined(vinecopulib_EXPORTS)
#        define  VINECOPULIB_EXPORT __declspec(dllexport)
#    else
#        define  VINECOPULIB_EXPORT __declspec(dllimport)
#    endif 
#else
#    define VINECOPULIB_EXPORT
#endif

// include all family headers not included yet (alphabetic order);
// this automatically includes the headers:
#include "bicop_archimedean.hpp"
#include "bicop_elliptical.hpp"
#include "bicop_parametric.hpp"
#include "bicop_kernel.hpp"
#include "bicop_bb1.hpp"
#include "bicop_bb6.hpp"
#include "bicop_bb7.hpp"
#include "bicop_bb8.hpp"
#include "bicop_clayton.hpp"
#include "bicop_frank.hpp"
#include "bicop_gaussian.hpp"
#include "bicop_gumbel.hpp"
#include "bicop_indep.hpp"
#include "bicop_joe.hpp"
#include "bicop_student.hpp"
#include "bicop_trafokernel.hpp"
