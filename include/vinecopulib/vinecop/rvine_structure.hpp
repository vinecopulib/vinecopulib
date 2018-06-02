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
    RVineStructure(const std::vector<size_t>& order,
                   const RVineMatrix<size_t>& struct_mat,
                   bool is_natural_order = false);

    std::vector<size_t> get_order() const;
    RVineMatrix<size_t> get_struct_matrix() const;
    RVineMatrix<size_t> get_max_matrix() const;
    RVineMatrix<size_t> get_needed_hfunc1() const;
    RVineMatrix<size_t> get_needed_hfunc2() const;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    size_t struct_matrix(size_t tree, size_t edge) const;
    size_t max_matrix(size_t tree, size_t edge) const;

protected:

    size_t find_trunc_lvl(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    std::vector<size_t> compute_order(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    RVineMatrix<size_t> compute_natural_order(const RVineMatrix<size_t>& mat) const;
    RVineMatrix<size_t> compute_struct_matrix(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    RVineMatrix<size_t> compute_dvine_struct_matrix() const;
    RVineMatrix<size_t> compute_max_matrix() const;
    RVineMatrix<size_t> compute_needed_hfunc1() const;
    RVineMatrix<size_t> compute_needed_hfunc2() const;


private:

    std::vector<size_t> order_;
    size_t d_;
    size_t trunc_lvl_;
    RVineMatrix<size_t> struct_mat_;
    RVineMatrix<size_t> max_mat_;
    RVineMatrix<size_t> needed_hfunc1_;
    RVineMatrix<size_t> needed_hfunc2_;
};

#include <vinecopulib/vinecop/implementation/rvine_structure.ipp>

}

