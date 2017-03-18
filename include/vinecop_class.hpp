// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <limits>

#include "bicop_class.hpp"
#include "rvine_matrix.hpp"

namespace vinecopulib
{
//! A class for vine copulas
    class Vinecop
    {
    public:
        Vinecop() {}
        Vinecop(int d);
        Vinecop(
                const std::vector<std::vector<BicopPtr>>& pair_copulas,
                const Eigen::MatrixXi& matrix
        );

        static std::vector<std::vector<BicopPtr>> make_pair_copula_store(int d);
        static Vinecop select(
                const Eigen::MatrixXd& data,
                std::vector<BicopFamily> family_set = bicop_families::all,
                std::string method = "mle",
                int truncation_level = std::numeric_limits<int>::max(),
                Eigen::MatrixXi matrix = Eigen::MatrixXi(0, 0),
                std::string selection_criterion = "bic",
                bool preselect_families = true,
                bool show_trace = false
        );

        BicopPtr get_pair_copula(int tree, int edge);
        BicopFamily get_family(int tree, int edge);
        std::vector<std::vector<BicopFamily>> get_all_families();
        int get_rotation(int tree, int edge);
        std::vector<std::vector<int>> get_all_rotations();
        std::vector<Eigen::MatrixXd> get_parameters(int tree, int edge);
        std::vector<std::vector<std::vector<Eigen::MatrixXd>>> get_all_parameters();
        Eigen::MatrixXi get_matrix() {return vine_matrix_.get_matrix();}

        Eigen::VectorXd pdf(const Eigen::MatrixXd& u);
        Eigen::MatrixXd simulate(int n);
        Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd& U);

    private:
        int d_;
        RVineMatrix vine_matrix_;
        std::vector<std::vector<BicopPtr>> pair_copulas_;
    };

    Eigen::VectorXi inverse_permutation(const Eigen::VectorXi& order);
    // reverse columns and rows of an Eigen::Matrix type object
    template<typename Mat>
    Mat to_upper_tri(Mat A) {return A.rowwise().reverse().colwise().reverse();}
}
