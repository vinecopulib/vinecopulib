// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

// interface specfifc #defines can be set here 
// (R package does: #define INTERFACED_FROM_R)


// interface specific headers
#ifdef INTERFACED_FROM_R
    #include <RcppThread.h>
#else
    #include <iostream>
#endif

// parallel backend
#ifdef INTERFACED_FROM_R
    namespace tools_parallel {
        typedef ThreadPool RcppThread::ThreadPool;
    }
#else
    #include <vinecopulib/misc/tools_parallel.hpp>
#endif

// for std::getenv() and std::atol()
#include <cstdlib>

namespace vinecopulib {

    namespace tools_interface {
        inline void print(std::string text)
        {
            #ifndef INTERFACED_FROM_R
                std::cout << text;
            #else
                RcppThread::Rcout << text;
            #endif
        }
    
        inline void check_user_interrupt(bool do_check = true)
        {
            if (do_check) {
                #ifdef INTERFACED_FROM_R
                    RcppThread::checkUserInterrupt();
                #endif
            }
        }
        
        inline size_t get_num_threads()
        {
            char* num_char = std::getenv("VINECOPULIB_NUM_THREADS");
            std::cout << "VINECOPULIB_NUM_THREADS = " << std::atoi(num_char) << "\n";
            // use at least one thread
            unsigned int num = std::max(std::atoi(num_char), 1);
            // don't use more threads than supported by the system
            num = std::min(num, std::thread::hardware_concurrency());
            
            return num;
        }
    }
}
