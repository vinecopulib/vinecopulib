// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <vinecopulib/vinecop/rvine_matrix.hpp>

namespace vinecopulib {

class RVineStructure {
public:
    RVineStructure() {}
    
    RVineStructure(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat,
        bool check = true);
    RVineStructure(const std::vector<size_t>& order,
                   bool check = true);
    RVineStructure(const std::vector<size_t>& order,
                   const size_t& trunc_lvl,
                   bool check = true);
    RVineStructure(const std::vector<size_t>& order,
                   const TriangularArray<size_t>& struct_mat,
                   bool is_natural_order = false,
                   bool check = true);

    size_t get_dim() const;
    size_t get_trunc_lvl() const;
    std::vector<size_t> get_order() const;
    TriangularArray<size_t> get_struct_matrix() const;
    TriangularArray<size_t> get_max_matrix() const;
    TriangularArray<size_t> get_needed_hfunc1() const;
    TriangularArray<size_t> get_needed_hfunc2() const;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    size_t struct_matrix(size_t tree, size_t edge) const;
    size_t max_matrix(size_t tree, size_t edge) const;

private:

    size_t find_trunc_lvl(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    std::vector<size_t> get_order(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    TriangularArray<size_t> to_rvine_matrix(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;

    TriangularArray<size_t> to_natural_order() const;
    TriangularArray<size_t> compute_dvine_struct_matrix() const;
    TriangularArray<size_t> compute_max_matrix() const;
    TriangularArray<size_t> compute_needed_hfunc1() const;
    TriangularArray<size_t> compute_needed_hfunc2() const;

    void check_if_quadratic(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_lower_tri(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_upper_tri() const;
    void check_columns() const;
    void check_antidiagonal() const;
    void check_proximity_condition() const;

    std::vector<size_t> order_;
    size_t d_;
    size_t trunc_lvl_;
    TriangularArray<size_t> struct_mat_;
    TriangularArray<size_t> max_mat_;
    TriangularArray<size_t> needed_hfunc1_;
    TriangularArray<size_t> needed_hfunc2_;
};

}

#include <vinecopulib/vinecop/implementation/rvine_structure.ipp>
