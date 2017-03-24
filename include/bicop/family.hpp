// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <vector>
#include <string>

//! @file family.hpp

namespace vinecopulib {
    /// identifies bivariate copula families.
    enum class BicopFamily
    {
        indep,         ///< Independence
        gaussian,      ///< Independence
        student,       ///< Independence
        clayton,       ///< Independence
        gumbel,        ///< Independence
        frank,         ///< Independence
        joe,           ///< Independence
        bb1,           ///< Independence
        bb6,           ///< Independence
        bb7,           ///< Independence
        bb8,           ///< Independence
        tll0           ///< Independence
    };

    std::string get_family_name(BicopFamily family);
    
    namespace bicop_families {
                
        const std::vector<BicopFamily> all = {
            BicopFamily::indep,
            BicopFamily::gaussian,
            BicopFamily::student,
            BicopFamily::clayton,
            BicopFamily::gumbel,
            BicopFamily::frank,
            BicopFamily::joe,
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8,
            BicopFamily::tll0
        };
        const std::vector<BicopFamily> parametric = {
            BicopFamily::indep,
            BicopFamily::gaussian,
            BicopFamily::student,
            BicopFamily::clayton,
            BicopFamily::gumbel,
            BicopFamily::frank,
            BicopFamily::joe,
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8
        };
        const std::vector<BicopFamily> nonparametric = {
            BicopFamily::indep,
            BicopFamily::tll0
        };
        const std::vector<BicopFamily> one_par = {
            BicopFamily::gaussian,
            BicopFamily::clayton,
            BicopFamily::gumbel,
            BicopFamily::frank,
            BicopFamily::joe,
        }; 
        const std::vector<BicopFamily> two_par = {
            BicopFamily::student,
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8
        };    
        const std::vector<BicopFamily> elliptical = {
            BicopFamily::gaussian, 
            BicopFamily::student
        };
        const std::vector<BicopFamily> archimedean = {
            BicopFamily::clayton, 
            BicopFamily::gumbel, 
            BicopFamily::frank, 
            BicopFamily::joe, 
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8
        };
        const std::vector<BicopFamily> BB = {
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8
        };        
        const std::vector<BicopFamily> rotationless = {
            BicopFamily::indep, 
            BicopFamily::gaussian, 
            BicopFamily::student, 
            BicopFamily::frank, 
            BicopFamily::tll0
        };
        const std::vector<BicopFamily> lt = {
            BicopFamily::clayton, 
            BicopFamily::bb1, 
            BicopFamily::bb7
        };  
        const std::vector<BicopFamily> ut = {
            BicopFamily::gumbel, 
            BicopFamily::joe, 
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8
        };      
        const std::vector<BicopFamily> itau = {
            BicopFamily::indep, 
            BicopFamily::gaussian, 
            BicopFamily::student, 
            BicopFamily::clayton, 
            BicopFamily::gumbel, 
            BicopFamily::frank, 
            BicopFamily::joe
        };
        
        const std::vector<BicopFamily> flip_by_rotation = {
            BicopFamily::clayton,
            BicopFamily::gumbel,
            BicopFamily::frank,
            BicopFamily::joe,
            BicopFamily::bb1, 
            BicopFamily::bb6, 
            BicopFamily::bb7, 
            BicopFamily::bb8
        };
        
    } // end of namespace BicopFamilies
    
} // end of namespace vinecopulib
