// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/factcop/class.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <stdexcept>
#include <vector>

namespace vinecopulib
{
    Factcop::Factcop(size_t d)
    {
        d_ = d;
        ngroups_ = 1;
        nfactors_ = 1;
        groups_ = Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Constant(d, 1);
        factors_ = Eigen::Matrix<bool, 1, 1>::Constant(true);

        // a single independence pair-copula
        pair_copulas_ = make_pair_copula_store(1, 1);
        pair_copulas_[0][0] = Bicop(BicopFamily::indep);
    }

    Factcop::Factcop(Eigen::Matrix<size_t, Eigen::Dynamic, 1> groups)
    {
        d_ = (size_t) groups.size();

        // get the number of unique groups
        std::vector<size_t> groups_stl(d_);
        Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&groups_stl[0], d_) = groups;
        ngroups_ = tools_stl::unique(groups_stl).size();

        nfactors_ = 1;
        groups_ = groups;
        factors_ = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(ngroups_,
                                                                    true);

        // all pair-copulas are independence
        pair_copulas_ = make_pair_copula_store(1, ngroups_);
        for (size_t group = 0; group < ngroups_; ++ group) {
            pair_copulas_[0][group] = Bicop(BicopFamily::indep);
        }

    }

    //! Initialize object for storing pair copulas
    //!
    //! @param f number of factors
    //! @param g number of groups
    //! @return A nested vector such that `pc_store[f][g]` contains a Bicop
    //!     object for the pair copula corresponding to factor `f` and group `g`.
    std::vector<std::vector<Bicop>> Factcop::make_pair_copula_store(
            size_t f,
            size_t g)
    {
        if (f < 1) {
            throw std::runtime_error("the number of factors should be larger than 0");
        }
        if (g < 1) {
            throw std::runtime_error("the number of groups should be larger than 0");
        }

        std::vector<std::vector<Bicop>> pc_store(f);
        for (size_t t = 0; t < f; ++t) {
            pc_store[t].resize(g);
        }

        return pc_store;
    }
}
