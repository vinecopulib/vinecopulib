// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/quickpool.hpp>
#include <mutex>
#include <condition_variable>

namespace vinecopulib {
namespace tools_thread {

// Derived ThreadPool class that extends quickpool::ThreadPool
class ThreadPool : public quickpool::ThreadPool
{
public:
    using quickpool::ThreadPool::ThreadPool; // Inherit constructors

    //! maps a function on a list of items, possibly running tasks in parallel.
    //! @param f Function to be mapped.
    //! @param items An object containing the items on which `f` shall be
    //!   mapped; must allow for `auto` loops (i.e., `std::begin(I)`/
    //!  `std::end(I)` must be defined).
    template <class F, class I>
    void map(F&& f, I&& items)
    {
        for (auto&& item : items) {
            this->push(f, item); // quickpool's push method
        }
    }
    //! waits for all jobs to finish.
    inline void join()
    {
        this->wait(); // Ensure all tasks are completed
        // The destructor will automatically stop the threads
    }
};

} // namespace tools_thread
} // namespace vinecopulib
