// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "misc/tools_stats.hpp"
#include "misc/tools_stl.hpp"
#include "misc/tools_c.h"

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
            x.col(j) = to_pseudo_obs_1d((Eigen::VectorXd) x.col(j), ties_method);

        return x;
    }

    Eigen::VectorXd to_pseudo_obs_1d(Eigen::VectorXd x, std::string ties_method)
    {
        int n = x.size();
        std::vector<double> xvec(x.data(), x.data() + n);
        auto order = tools_stl::get_order(xvec);
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

    double pairwise_ktau(Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        double tau;
        int n = u.rows();
        int two = 2;
        ktau_matrix(u.data(), &two, &n, &tau);
        return tau;
    }

    double pairwise_cor(const Eigen::Matrix<double, Eigen::Dynamic, 2>& z)
    {
        double rho;
        auto x = z.rowwise() - z.colwise().mean();
        Eigen::MatrixXd sigma = x.adjoint() * x;
        rho = sigma(1,0) / sqrt(sigma(0,0) * sigma(1,1));

        return rho;
    }

    double pairwise_hoeffd(Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        int n = x.rows();

        // Compute the ranks
        auto R = to_pseudo_obs(x);
        R = (n+1.0)*R;

        // Compute Q, with Qi the number of points with both columns less than
        // their ith value
        Eigen::VectorXd Q(n);
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = Eigen::MatrixXd::Ones(n, 2);
        for(int i=0; i<n; i++) {
            tmp.col(0) = Eigen::VectorXd::Constant(n,x(i,0));
            tmp.col(1) = Eigen::VectorXd::Constant(n,x(i,1));
            tmp = (x-tmp).unaryExpr([](double v){
                double res = 0.0;
                if(v < 0.0) {
                    res = 1.0;
                }
                return res;
            });
            Q(i) = tmp.rowwise().prod().sum();
        }

        Eigen::Matrix<double, Eigen::Dynamic, 2> ones = Eigen::MatrixXd::Ones(n, 2);
        double A = (R-ones).cwiseProduct(R-2*ones).rowwise().prod().sum();
        double B = (R-2*ones).rowwise().prod().cwiseProduct(Q).sum();
        double C = Q.cwiseProduct(Q-ones.col(0)).sum();

        double D = (A - 2*(n-2)*B + (n-2)*(n-3)*C);
        D /= (n*(n-1)*(n-2)*(n-3)*(n-4));

        return D;
    }
}
