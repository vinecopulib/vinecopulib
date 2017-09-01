// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/class.hpp>

namespace vinecopulib
{
    //! @brief A class for factor copula models
    //!
    class Factcop
    {
    public:
        // Constructors
        Factcop() {}
        Factcop(size_t d);
        Factcop(Eigen::Matrix<size_t, Eigen::Dynamic, 1> groups);

        std::vector<std::vector<Bicop>> make_pair_copula_store(size_t f,
                                                               size_t g);

    private:
        size_t d_;
        size_t ngroups_;
        size_t nfactors_;
        Eigen::Matrix<size_t, Eigen::Dynamic, 1> groups_;
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> factors_;
        std::vector<std::vector<Bicop>> pair_copulas_;
    };

}
