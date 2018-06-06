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
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat);
    RVineStructure(const std::vector<size_t>& order);
    RVineStructure(const std::vector<size_t>& order, const size_t& trunc_lvl);
    RVineStructure(const std::vector<size_t>& order,
                   const RVineMatrix<size_t>& struct_mat,
                   bool is_natural_order = false);

    size_t get_dim() const;
    size_t get_trunc_lvl() const;
    std::vector<size_t> get_order() const;
    RVineMatrix<size_t> get_struct_matrix() const;
    RVineMatrix<size_t> get_max_matrix() const;
    RVineMatrix<size_t> get_needed_hfunc1() const;
    RVineMatrix<size_t> get_needed_hfunc2() const;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    size_t struct_matrix(size_t tree, size_t edge) const;
    size_t max_matrix(size_t tree, size_t edge) const;

private:

    size_t find_trunc_lvl(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    std::vector<size_t> get_order(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;

    template<class T> RVineMatrix<size_t> to_natural_order(const T& mat) const;
    RVineMatrix<size_t> compute_dvine_struct_matrix() const;
    RVineMatrix<size_t> compute_max_matrix() const;
    RVineMatrix<size_t> compute_needed_hfunc1() const;
    RVineMatrix<size_t> compute_needed_hfunc2() const;

    void check_if_quadratic(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_lower_tri(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_upper_tri(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_antidiagonal(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;

    std::vector<size_t> order_;
    size_t d_;
    size_t trunc_lvl_;
    RVineMatrix<size_t> struct_mat_;
    RVineMatrix<size_t> max_mat_;
    RVineMatrix<size_t> needed_hfunc1_;
    RVineMatrix<size_t> needed_hfunc2_;
};

template<class T> RVineMatrix<size_t> RVineStructure::to_natural_order(const T& mat) const
{
    // create vector of new variable labels
    auto order = tools_stl::get_order(get_order());

    // copy upper triangle and relabel to natural order
    RVineMatrix<size_t> struct_mat(d_, trunc_lvl_);
    for (size_t j = 0; j < d_ - 1; j++) {
        for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
            struct_mat(i, j) = order[mat(i, j) - 1] + 1;
        }
    }

    return struct_mat;
}

}

#include <vinecopulib/vinecop/implementation/rvine_structure.ipp>
