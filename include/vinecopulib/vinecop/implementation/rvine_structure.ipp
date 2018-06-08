// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

inline RVineStructure::RVineStructure(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat,
    bool check)
{
    d_ = mat.cols();
    if (check) {
        check_if_quadratic(mat);
        check_lower_tri(mat);
    }

    order_ = get_order(mat);
    if (check)
        check_antidiagonal();

    trunc_lvl_ = find_trunc_lvl(mat);
    struct_mat_ = to_rvine_matrix(mat);
    if (check)
        check_upper_tri();

    struct_mat_ = to_natural_order();
    if (check)
        check_columns();

    max_mat_ = compute_max_matrix();
    if (check)
        check_proximity_condition();

    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
}

inline RVineStructure::RVineStructure(
    const std::vector<size_t>& order,
    bool check) : RVineStructure(order, order.size() - 1, check) {}

inline RVineStructure::RVineStructure(
    const std::vector<size_t>& order, const size_t& trunc_lvl, bool check)
{
    d_ = order.size();
    order_ = order;

    if (check)
        check_antidiagonal();

    if (trunc_lvl > 0) {
        if (trunc_lvl > d_ - 1) {
            trunc_lvl_ = d_ - 1;
        } else {
            trunc_lvl_ = trunc_lvl;
        }

        struct_mat_ = compute_dvine_struct_matrix();
        max_mat_ = compute_max_matrix();
        needed_hfunc1_ = compute_needed_hfunc1();
        needed_hfunc2_ = compute_needed_hfunc2();
    }
}

inline RVineStructure::RVineStructure(
    const std::vector<size_t>& order,
    const TriangularArray<size_t>& struct_mat,
    bool is_natural_order,
    bool check)
{
    d_ = order.size();
    if (check & (struct_mat.get_dim() != d_)) {
        throw std::runtime_error("order and struct_mat have "
                                     "incompatible dimensions");
    }

    order_ = order;

    if (check)
        check_antidiagonal();

    trunc_lvl_ = struct_mat.get_trunc_lvl();
    if (trunc_lvl_ > 0) {

        struct_mat_ = struct_mat;
        if (check)
            check_upper_tri();

        if (!is_natural_order)
            struct_mat_ = to_natural_order();
        if (check)
            check_columns();

        max_mat_ = compute_max_matrix();
        if (check)
            check_proximity_condition();

        needed_hfunc1_ = compute_needed_hfunc1();
        needed_hfunc2_ = compute_needed_hfunc2();
    }
}

inline size_t RVineStructure::get_dim() const
{
    return d_;
}

inline size_t RVineStructure::get_trunc_lvl() const
{
    return trunc_lvl_;
}

inline std::vector<size_t> RVineStructure::get_order() const 
{
    return order_;
}

inline TriangularArray<size_t> RVineStructure::get_struct_matrix() const 
{
    return struct_mat_;
}

inline TriangularArray<size_t> RVineStructure::get_max_matrix() const 
{
    return max_mat_;
}

inline TriangularArray<size_t> RVineStructure::get_needed_hfunc1() const 
{
    return needed_hfunc1_;
}

inline TriangularArray<size_t> RVineStructure::get_needed_hfunc2() const
{
    return needed_hfunc2_;
}

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
    for (size_t i = 0; i < trunc_lvl_; ++i) {
        for (size_t j = 0; j < d_ - i - 1; ++j) {
            matrix(i, j) = order_[struct_mat_(i, j) - 1];
        }
    }
    for (size_t i = 0; i < d_; ++i) {
        matrix(d_ - i - 1, i) = order_[d_ - i - 1];
    }
    return matrix;
}

inline size_t RVineStructure::find_trunc_lvl(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    size_t trunc_lvl;
    size_t d = mat.cols();

    std::stringstream problem;
    problem << "not a valid R-vine matrix: " <<
            "a row with a 0 above the diagonal contains one or more " <<
            "non-zero values.";

    for (trunc_lvl = d - 1; trunc_lvl > 0; trunc_lvl--) {
        std::vector<size_t> row_vec(d - trunc_lvl);
        Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&row_vec[0],
                                                      d - trunc_lvl) =
            mat.row(trunc_lvl - 1).head(d - trunc_lvl);

        if (*(std::min_element(row_vec.begin(), row_vec.end())) != 0)
            break;
    }

    return trunc_lvl;
}

inline std::vector<size_t> RVineStructure::get_order(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    std::vector<size_t> order(d_);
    for (size_t i = 0; i < d_; i++)
        order[i] = mat(i, d_ - i - 1);

    return order;
}

inline TriangularArray<size_t> RVineStructure::to_rvine_matrix(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    // copy upper triangle
    TriangularArray<size_t> struct_mat(d_, trunc_lvl_);
    for (size_t j = 0; j < d_ - 1; j++) {
        for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
            struct_mat(i, j) = mat(i, j);
        }
    }

    return struct_mat;
}

inline TriangularArray<size_t> RVineStructure::to_natural_order() const
{
    // create vector of new variable labels
    auto order = tools_stl::get_order(get_order());

    // relabel to natural order
    TriangularArray<size_t> struct_mat(d_, trunc_lvl_);
    for (size_t j = 0; j < d_ - 1; j++) {
        for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
            struct_mat(i, j) = order[struct_mat_(i, j) - 1] + 1;
        }
    }

    return struct_mat;
}

inline TriangularArray<size_t> RVineStructure::compute_dvine_struct_matrix() const
{
    TriangularArray<size_t> struct_mat(d_, trunc_lvl_);
    for (size_t j = 0; j < d_ - 1; j++) {
        for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
            struct_mat(i, j) = d_ - i - j - 1;
        }
    }

    return struct_mat;
}

inline TriangularArray<size_t> RVineStructure::compute_max_matrix() const
{
    TriangularArray<size_t> max_mat = struct_mat_;
    for (size_t j = 0; j < d_ - 1; j++) {
        for (size_t i = 1; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
            max_mat(i, j) = std::max(struct_mat_(i, j), max_mat(i - 1, j));
        }
    }

    return max_mat;
}

inline TriangularArray<size_t> RVineStructure::compute_needed_hfunc1() const
{
    TriangularArray<size_t> needed_hfunc1(d_, trunc_lvl_);

    for (size_t i = 0; i < std::min(d_ - 2, trunc_lvl_ - 1); i++) {
        for (size_t j = 0; j < d_ - 2 - i; j++) {
            if (struct_mat_(i + 1, j) != max_mat_(i + 1, j))
                needed_hfunc1(i, d_ - max_mat_(i + 1, j)) = 1;
        }
    }

    return needed_hfunc1;
}

inline TriangularArray<size_t> RVineStructure::compute_needed_hfunc2() const
{
    TriangularArray<size_t> needed_hfunc2(d_, trunc_lvl_);

    for (size_t i = 0; i < std::min(d_ - 2, trunc_lvl_ - 1); i++) {
        for (size_t j = 0; j < d_ - 2 - i; j++) {
            needed_hfunc2(i, j) = 1;
            if (struct_mat_(i + 1, j) == max_mat_(i + 1, j))
                needed_hfunc2(i, d_ - max_mat_(i + 1, j)) = 1;
        }
    }

    return needed_hfunc2;
}

inline void RVineStructure::check_if_quadratic(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    std::string problem = "must be quadratic.";
    if (mat.rows() != mat.cols()) {
        throw std::runtime_error("not a valid R-vine matrix: " + problem);
    }
}

inline void RVineStructure::check_lower_tri(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
    std::string problem = "the lower right triangle must only contain zeros";
    size_t sum_lwr = 0;
    for (size_t j = 1; j < d_; ++j) {
        sum_lwr += mat.block(d_ - j, j, j, 1).array().sum();
        if (sum_lwr != 0) {
            throw std::runtime_error("not a valid R-vine matrix: " + problem);
        }
    }
}

inline void RVineStructure::check_upper_tri() const
{
    std::string problem;
    problem += "the upper left triangle can only contain numbers ";
    problem += "between 1 and d (number of variables).";

    for (size_t j = 0; j < d_ - 1; ++j) {
        auto col_vec = struct_mat_[j];
        auto minmax_in_col = std::minmax_element(col_vec.begin(),
                                                 col_vec.end());
        if ((*(minmax_in_col.first) < 1) | (*(minmax_in_col.second) > d_)) {
            throw std::runtime_error("not a valid R-vine matrix: " + problem);
        }
    }
}

inline void RVineStructure::check_columns() const
{
    std::string problem;
    problem += "the antidiagonal entry of a column must not be ";
    problem += "contained in any column further to the right; ";
    problem += "the entries of a column must be contained ";
    problem += "in all columns to the left.";

    // check that column j only contains unique indices in 1:(d - j).
    for (size_t j = 0; j < d_ - 1; ++j) {
        auto col_vec = struct_mat_[j];
        std::sort(col_vec.begin(), col_vec.end());
        size_t unique_in_col = std::unique(col_vec.begin(), col_vec.end())
                               - col_vec.begin();
        if ((!tools_stl::is_member(col_vec, tools_stl::seq_int(1, d_ - j))) |
            (unique_in_col != col_vec.size())) {
            throw std::runtime_error("not a valid R-vine matrix: " + problem);
        }
    }
}

inline void RVineStructure::check_antidiagonal() const
{
    std::string problem;
    problem += "the order/antidiagonal must contain the numbers ";
    problem += "1, ..., d (the number of variables)";
    if (!tools_stl::is_same_set(order_, tools_stl::seq_int(1, d_))) {
        throw std::runtime_error("not a valid R-vine matrix: " + problem);
    }
}


inline void RVineStructure::check_proximity_condition() const
{
    for (size_t t = 1; t < trunc_lvl_; ++t) {
        for (size_t e = 0; e < d_ - t - 1; ++e) {
            std::vector <size_t> target_set(t + 1), test_set(t + 1);
            // conditioning set
            for (size_t i = 0; i < t; i++) {
                target_set[i] = struct_mat_(i, e);
                test_set[i] = struct_mat_(i, d_ - max_mat_(t, e));
            }
            // non-diagonal conditioned variable
            target_set[t] = struct_mat_(t, e);
            // diagonal conditioned variable in other column
            test_set[t] = max_mat_(t, e);

            if (!tools_stl::is_same_set(target_set, test_set)) {
                std::stringstream problem;
                problem << "not a valid R-vine matrix: " <<
                        "proximity condition violated; " <<
                        "cannot extract conditional distribution (" <<
                        target_set[t] << " | ";
                for (size_t i = 0; i < t - 1; ++i) {
                    problem << order_[target_set[i] - 1] << ", ";
                }
                problem << order_[target_set[t - 1] - 1]
                        << ") from pair-copulas.";
                throw std::runtime_error(problem.str().c_str());
            }
        }
    }
}

}
