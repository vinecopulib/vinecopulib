// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <limits>

#include "bicop/class.hpp"
#include "vinecop/rvine_matrix.hpp"

namespace vinecopulib
{
//! A class for vine copulas
    class Vinecop
    {
    public:
        Vinecop() {}
        Vinecop(int d);
        Vinecop(
                const std::vector<std::vector<Bicop>>& pair_copulas,
                const Eigen::MatrixXi& matrix
        );

        static std::vector<std::vector<Bicop>> make_pair_copula_store(int d);
        static Vinecop select(
                const Eigen::MatrixXd& data,
                std::vector<BicopFamily> family_set = bicop_families::all,
                std::string method = "mle",
                int truncation_level = std::numeric_limits<int>::max(),
                double threshold = 0.0,
                Eigen::MatrixXi matrix = Eigen::MatrixXi(0, 0),
                std::string tree_criterion = "tau",
                std::string selection_criterion = "bic",
                bool preselect_families = true,
                bool show_trace = false
        );

        Bicop get_pair_copula(int tree, int edge) const;
        BicopFamily get_family(int tree, int edge) const;
        std::vector<std::vector<BicopFamily>> get_all_families() const;
        int get_rotation(int tree, int edge) const;
        std::vector<std::vector<int>> get_all_rotations() const;
        Eigen::VectorXd get_parameters(int tree, int edge) const;
        std::vector<std::vector<Eigen::VectorXd>> get_all_parameters() const;
        Eigen::MatrixXi get_matrix() const;

        Eigen::VectorXd pdf(const Eigen::MatrixXd& u);
        Eigen::MatrixXd simulate(int n);
        Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd& U);

    private:
        int d_;
        RVineMatrix vine_matrix_;
        std::vector<std::vector<Bicop>> pair_copulas_;
    };

    Eigen::VectorXi inverse_permutation(const Eigen::VectorXi& order);
    // reverse columns and rows of an Eigen::Matrix type object
    template<typename Mat>
    Mat to_upper_tri(Mat A) {return A.rowwise().reverse().colwise().reverse();}
}
