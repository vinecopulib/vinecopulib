// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "tools_stats.hpp"
#include <exception>
#include <algorithm>

namespace tools_stats {
    
    //! Empirical probability integral transform
    //!
    //! Gives pseudo-observations from the copula by applying the empirical
    //! distribution function (scaled by n + 1) to each margin/column.
    //!
    //! @param x a matrix of real numbers.
    //! @param ties_method indicates how to treat ties; same as in R, see
    //! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
    //! @return Psuedo-observations of the copula, i.e. F_X(X) (column-wise)
    Eigen::MatrixXd to_pseudo_obs(Eigen::MatrixXd x, std::string ties_method)
    {
        for (unsigned int j = 0; j < x.cols(); ++j)
            x.col(j) = to_pseudo_obs((Eigen::VectorXd) x.col(j), ties_method);

        return x;
    }

    Eigen::VectorXd to_pseudo_obs(Eigen::VectorXd x, std::string ties_method)
    {
        int n = x.size();
        std::vector<double> xvec(x.data(), x.data() + n);
        auto order = get_order(xvec);
        if (ties_method == "first") {
            for (auto i : order)
                x[order[i]] = i + 1;
        } else if (ties_method == "average") {
            for (int i = 0, reps; i < n; i += reps) {
                // find replications
                reps = 1;
                while ((i + reps < n) && (x[order[i]] == x[order[i + reps]]))
                    ++reps;
                // assign average rank of the tied values
                for (int k = 0; k < reps; ++k)
                    x[order[i + k]] = i + 1 + (reps - 1) / 2.0;
            }
        } else if (ties_method == "random") {
            // set up random number generator
            std::random_device rd;
            std::default_random_engine gen(rd());
            auto sim = [&] (int m) {
                std::uniform_int_distribution<> distr(0, m - 1);
                return distr(gen);
            };
            for (int i = 0, reps; i < n; i += reps) {
                // find replications
                reps = 1;
                while ((i + reps < n) && (x[order[i]] == x[order[i + reps]]))
                    ++reps;
                // assign random rank between ties
                std::vector<int> rvals(reps);
                std::iota(rvals.begin(), rvals.end(), 0);  // 0, 1, 2, ...
                std::random_shuffle(rvals.begin(), rvals.end(), sim);
                for (int k = 0; k < reps; ++k)
                    x[order[i + k]] = i + 1 + rvals[k];
            }
        } else {
            std::stringstream msg;
            msg << "unknown ties method (" << ties_method << ")";
            throw std::runtime_error(msg.str().c_str());
        }

        return x / (x.size() + 1.0);
    }
}
