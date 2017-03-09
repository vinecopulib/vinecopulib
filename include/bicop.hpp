/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see "http://www.gnu.org/licenses/".
*/

#ifndef VINECOPULIB_BICOP_FAMILIES_HPP
#define VINECOPULIB_BICOP_FAMILIES_HPP

    // include all family headers not included yet (alphabetic order);
    // this automatically includes the headers:
    //   bicop_class
    //   bicop_archimedean
    //   bicop_elliptical
    //   bicop_parametric
    //   bicop_kernel
    //   bicop_interpolation
    #include "bicop_bb1.hpp"
    #include "bicop_bb6.hpp"
    #include "bicop_clayton.hpp"
    #include "bicop_frank.hpp"
    #include "bicop_gauss.hpp"
    #include "bicop_gumbel.hpp"
    #include "bicop_indep.hpp"
    #include "bicop_joe.hpp"
    #include "bicop_student.hpp"
    #include "bicop_trafokernel.hpp"

#endif
