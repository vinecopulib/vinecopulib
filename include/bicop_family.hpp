// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <vector>

namespace vinecopulib {
    //! Bivariate copula families
    enum class BicopFamily
    {
        Indep,
        Gaussian,
        Student,
        Clayton,
        Gumbel,
        Frank,
        Joe,
        BB1,
        BB6,
        BB7,
        BB8,
        TLL0
    };
    
    namespace bicop_families {
                
        const std::vector<BicopFamily> all = {
            BicopFamily::Indep,
            BicopFamily::Gaussian,
            BicopFamily::Student,
            BicopFamily::Clayton,
            BicopFamily::Gumbel,
            BicopFamily::Frank,
            BicopFamily::Joe,
            BicopFamily::BB1, 
            BicopFamily::BB6, 
            BicopFamily::BB7, 
            BicopFamily::BB8,
            BicopFamily::TLL0
        };
        const std::vector<BicopFamily> parametric = {
            BicopFamily::Gaussian,
            BicopFamily::Student,
            BicopFamily::Clayton,
            BicopFamily::Gumbel,
            BicopFamily::Frank,
            BicopFamily::Joe,
            BicopFamily::BB1, 
            BicopFamily::BB6, 
            BicopFamily::BB7, 
            BicopFamily::BB8
        };
        const std::vector<BicopFamily> nonparametric = {
            BicopFamily::Indep,
            BicopFamily::TLL0
        };
        const std::vector<BicopFamily> one_par = {
            BicopFamily::Gaussian,
            BicopFamily::Clayton,
            BicopFamily::Gumbel,
            BicopFamily::Frank,
            BicopFamily::Joe,
        }; 
        const std::vector<BicopFamily> two_par = {
            BicopFamily::Student,
            BicopFamily::BB1, 
            BicopFamily::BB6, 
            BicopFamily::BB7, 
            BicopFamily::BB8
        };    
        const std::vector<BicopFamily> elliptical = {
            BicopFamily::Gaussian, 
            BicopFamily::Student
        };
        const std::vector<BicopFamily> archimedean = {
            BicopFamily::Clayton, 
            BicopFamily::Gumbel, 
            BicopFamily::Frank, 
            BicopFamily::Joe, 
            BicopFamily::BB1, 
            BicopFamily::BB6, 
            BicopFamily::BB7, 
            BicopFamily::BB8
        };
        const std::vector<BicopFamily> BB = {
            BicopFamily::BB1, 
            BicopFamily::BB6, 
            BicopFamily::BB7, 
            BicopFamily::BB8
        };        
        const std::vector<BicopFamily> rotationless = {
            BicopFamily::Indep, 
            BicopFamily::Gaussian, 
            BicopFamily::Student, 
            BicopFamily::Frank, 
            BicopFamily::TLL0
        };
        const std::vector<BicopFamily> lt = {
            BicopFamily::Clayton, 
            BicopFamily::BB1, 
            BicopFamily::BB7
        };  
        const std::vector<BicopFamily> ut = {
            BicopFamily::Gumbel, 
            BicopFamily::Joe, 
            BicopFamily::BB1, 
            BicopFamily::BB6, 
            BicopFamily::BB7, 
            BicopFamily::BB8
        };      
        const std::vector<BicopFamily> itau = {
            BicopFamily::Indep, 
            BicopFamily::Gaussian, 
            BicopFamily::Student, 
            BicopFamily::Clayton, 
            BicopFamily::Gumbel, 
            BicopFamily::Frank, 
            BicopFamily::Joe
        };
        
    } // end of namespace BicopFamilies
    
} // end of namespace vinecopulib
