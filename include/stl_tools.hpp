// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

namespace stl_tools {

    template<typename T> bool is_member(T element, std::vector<T> set) {
        return std::find(set.begin(), set.end(), element) != set.end(); }

    template<class T>
    std::vector<T> intersect(std::vector<T> x, std::vector<T> y)
    {
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
        std::vector<T> common;
        std::set_intersection(
                x.begin(), x.end(),
                y.begin(), y.end(),
                std::back_inserter(common)
        );

        return common;
    }

    template<class T>
    T find_position(T x, std::vector<T> vec)
    {
        return std::distance(vec.begin(), std::find(vec.begin(), vec.end(), x));
    }

    template<class T>
    std::vector<T> set_diff(std::vector<T> x, std::vector<T> y)
    {
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
        std::vector<T> different;
        std::set_difference(
                x.begin(), x.end(),
                y.begin(), y.end(),
                std::back_inserter(different)
        );

        return different;
    }

    template<class T>
    std::vector<T> cat(std::vector<T> x, const std::vector<T>& y)
    {
        x.reserve(x.size() + y.size());
        x.insert(x.end(), y.begin(), y.end());
        return x;
    }

    template<class T>
    std::vector<T> cat(T x, const std::vector<T>& y)
    {
        std::vector<T> out(1);
        out[0] = x;
        out.reserve(1 + y.size());
        out.insert(out.end(), y.begin(), y.end());
        return out;
    }

    template<class T>
    void reverse(std::vector<T>& x)
    {
        std::reverse(x.begin(), x.end());
    }

    template<class T>
    bool is_same_set(std::vector<T> x, std::vector<T> y)
    {
        auto z = intersect(x, y);
        return ((z.size() == x.size()) & (z.size() == y.size()));
    }

    //! Integer sequence starting at 1
    inline std::vector<int> seq_int(int from, int length)
    {
        std::vector<int> seq(length);
        std::iota(seq.begin(), seq.end(), from);
        return seq;
    }
}
