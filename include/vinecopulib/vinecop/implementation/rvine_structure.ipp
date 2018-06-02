// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_parallel.hpp>

namespace vinecopulib {

inline RVineStructure::RVineStructure(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
    d_ = mat.cols();

    trunc_lvl_ = find_trunc_lvl(mat);
    order_ = compute_order(mat);
    struct_mat_ = compute_struct_matrix(mat);
    max_mat_ = compute_max_matrix();
    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
}

inline RVineStructure::RVineStructure(const std::vector<size_t>& order)
{

    d_ = order.size();
    order_ = order;
    trunc_lvl_ = d_;
    struct_mat_ = compute_dvine_struct_matrix();
    max_mat_ = compute_max_matrix();
    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();

}

inline RVineStructure::RVineStructure(
    const std::vector<size_t>& order,
    const RVineMatrix<size_t>& struct_mat,
    bool is_natural_order)
{
    d_ = order.size();
    order_ = order;
    trunc_lvl_ = d_;
    if (is_natural_order) {
        struct_mat_ = struct_mat;
    } else {
        struct_mat_ = compute_natural_order(struct_mat);
    }
    max_mat_ = compute_max_matrix();
    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
}

inline std::vector<size_t> RVineStructure::get_order() const 
{
    return order_;
}

inline RVineMatrix<size_t> RVineStructure::get_struct_matrix() const 
{
    return struct_mat_;
}

inline RVineMatrix<size_t> RVineStructure::get_max_matrix() const 
{
    return max_mat_;
}

inline RVineMatrix<size_t> RVineStructure::get_needed_hfunc1() const 
{
    return needed_hfunc1_;
}

inline RVineMatrix<size_t> RVineStructure::get_needed_hfunc2() const {return needed_hfunc2_;}

inline size_t RVineStructure::struct_matrix(size_t tree, size_t edge) const 
{
    return struct_mat_(tree, edge);
}

inline size_t RVineStructure::max_matrix(size_t tree, size_t edge) const {
    return max_mat_(tree, edge);
}

inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> RVineStructure::get_matrix() const 
{
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix(d_, d_);
    matrix.fill(0);
    for (size_t i = 0; i < d_; ++i) {
        for (size_t j = 0; j < d_ - i - 1; ++j) {
            matrix(i, j) = order_[struct_mat_(i, j) - 1];
        }
        matrix(d_ - i - 1, i) = order_[d_ - i - 1];
    }
    return matrix;
}

inline size_t RVineStructure::find_trunc_lvl(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    size_t trunc_lvl; ;
    for (trunc_lvl = mat.cols() - 1; trunc_lvl > 0; trunc_lvl--) {
        if (mat(trunc_lvl - 1, 0) != 0)
            break;
    }

    return trunc_lvl;
}

inline std::vector<size_t> RVineStructure::compute_order(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    std::vector<size_t> order(d_);
    for (size_t i = 0; i < d_; i++)
        order[i] = mat(i, d_ - i - 1);

    return order;
}

inline RVineMatrix<size_t> RVineStructure::compute_natural_order(
    const RVineMatrix<size_t>& mat) const
{
    // create vector of new variable labels
    auto order = tools_stl::get_order(get_order());

    // copy upper triangle and relabel to natural order
    RVineMatrix<size_t> struct_mat(d_);
    for (size_t i = 0; i < d_ - 1; i++) {
        for (size_t j = 0; j < d_ - 1 - i; j++) {
            struct_mat(i, j) = order[mat(i, j) - 1] + 1;
        }
    }

    return struct_mat;
}

inline RVineMatrix<size_t> RVineStructure::compute_struct_matrix(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    // create vector of new variable labels
    auto order = tools_stl::get_order(get_order());

    // copy upper triangle and relabel to natural order
    RVineMatrix<size_t> struct_mat(d_);
    for (size_t i = 0; i < d_ - 1; i++) {
        for (size_t j = 0; j < d_ - 1 - i; j++) {
            struct_mat(i, j) = order[mat(i, j) - 1] + 1;
        }
    }

    return struct_mat;
}

inline RVineMatrix<size_t> RVineStructure::compute_dvine_struct_matrix() const
{
    RVineMatrix<size_t> struct_mat(d_);
    for (size_t i = 0; i < d_ - 1; i++) {
        for (size_t j = 0; j < d_ - 1 - i; j++) {
            struct_mat(i, j) = d_ - i - j - 1;
        }
    }

    return struct_mat;
}

inline RVineMatrix<size_t> RVineStructure::compute_max_matrix() const
{
    RVineMatrix<size_t> max_mat = struct_mat_;
    for (size_t j = 0; j < d_ - 1; j++) {
        for (size_t i = 1; i < d_ - 1 - j; i++) {
            max_mat(i, j) = std::max(struct_mat_(i, j), max_mat(i - 1, j));
        }
    }

    return max_mat;
}

inline RVineMatrix<size_t> RVineStructure::compute_needed_hfunc1() const
{
    RVineMatrix<size_t> needed_hfunc1(d_);

    for (size_t i = 0; i < d_ - 2; i++) {
        for (size_t j = 0; j < d_ - 2 - i; j++) {
            if (struct_mat_(i + 1, j) != max_mat_(i + 1, j))
                needed_hfunc1(i, d_ - max_mat_(i + 1, j)) = 1;
        }
    }

    return needed_hfunc1;
}

inline RVineMatrix<size_t> RVineStructure::compute_needed_hfunc2() const
{
    RVineMatrix<size_t> needed_hfunc2(d_);

    for (size_t i = 0; i < d_ - 2; i++) {
        for (size_t j = 0; j < d_ - 2 - i; j++) {
            needed_hfunc2(i, j) = 1;
            if (struct_mat_(i + 1, j) == max_mat_(i + 1, j))
                needed_hfunc2(i, d_ - max_mat_(i + 1, j)) = 1;
        }
    }

    return needed_hfunc2;
}

}
