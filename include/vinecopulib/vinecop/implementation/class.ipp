// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_interface.hpp>

#include <stdexcept>

namespace vinecopulib {

//! creates a D-vine on `d` variables with all pair-copulas set to
//! independence.
//! @param d the dimension (= number of variables) of the model.
inline Vinecop::Vinecop(const size_t d)
{
    d_ = d;

    // 0-truncated D-vine with variable order (1, ..., d)
    vine_struct_ = RVineStructure(tools_stl::seq_int(1, d),
                                  static_cast<size_t>(0));

    // pair_copulas_ empty = everything independence 
    threshold_ = 0.0;
    loglik_ = NAN;
}

//! creates a vine copula with structure specified by an RVineStructure object;
//! all pair-copulas are set to independence.
//! @param vine_struct an RVineStructure object representing the structure of
//! the vine.
inline Vinecop::Vinecop(const RVineStructure &vine_struct)
{
    d_ = vine_struct.get_dim();
    vine_struct_ = vine_struct;
    // pair_copulas_ empty = everything independence
    threshold_ = 0.0;
    loglik_ = NAN;
}

//! creates a vine copula with structure specified by an R-vine matrix; all
//! pair-copulas are set to independence.
//! @param matrix an R-vine matrix.
//! @param check_matrix whether to check if `matrix` is a valid R-vine
//!     matrix.
inline Vinecop::Vinecop(
    const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
    const bool check_matrix) :
    Vinecop(RVineStructure(matrix, check_matrix)) {}

//! creates a vine copula with structure specified by an R-vine matrix; all
//! pair-copulas are set to independence.
//! @param order the order of the variables in the vine structure, see
//! RVineStructure's corresponding constructor.
//! @param struct_array a triangular array object specifying the vine structure,
//! see RVineStructure's corresponding constructor.
//! @param check_array whether `order` and `struct_array` shall be checked
//! for validity.
inline Vinecop::Vinecop(const std::vector<size_t> &order,
        const TriangularArray<size_t> &struct_array,
        const bool check_array)  :
    Vinecop(RVineStructure(order, struct_array, false, check_array)) {}

//! creates an arbitrary vine copula model.
//! @param pair_copulas Bicop objects specifying the pair-copulas, see
//!     make_pair_copula_store().
//! @param vine_struct an RVineStructure object specifying the vine structure.
inline Vinecop::Vinecop(const std::vector<std::vector<Bicop>> &pair_copulas,
                        const RVineStructure &vine_struct)
{

    d_ = vine_struct.get_dim();
    vine_struct_ = vine_struct;

    check_pair_copulas_rvine_structure(pair_copulas);

    pair_copulas_ = pair_copulas;
    threshold_ = 0.0;
    loglik_ = NAN;
}

//! creates an arbitrary vine copula model.
//! @param pair_copulas Bicop objects specifying the pair-copulas, see
//!     make_pair_copula_store().
//! @param matrix an R-vine matrix specifying the vine structure.
//! @param check_matrix whether to check if `matrix` is a valid R-vine
//!     matrix.
inline Vinecop::Vinecop(const std::vector<std::vector<Bicop>> &pair_copulas,
                        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
                        const bool check_matrix) :
    Vinecop(pair_copulas,
            RVineStructure(matrix, check_matrix)) {}

//! creates an arbitrary vine copula model.
//! @param pair_copulas Bicop objects specifying the pair-copulas, see
//!     make_pair_copula_store().
//! @param order the order of the variables in the vine structure, see
//! RVineStructure's corresponding constructor.
//! @param struct_array a triangular array object specifying the vine structure,
//! see the corresponding constructor of RVineStructure.
//! @param check_array whether `order` and `struct_array` shall be checked
//! for validity.
inline Vinecop::Vinecop(const std::vector<std::vector<Bicop>> &pair_copulas,
                        const std::vector<size_t> &order,
                        const TriangularArray<size_t> &struct_array,
                        const bool check_array) :
    Vinecop(pair_copulas,
            RVineStructure(order, struct_array, false, check_array)) {}

//! creates from a boost::property_tree::ptree object
//! @param input the boost::property_tree::ptree object to convert from
//! (see to_ptree() for the structure of the input).
//! @param check_matrix whether to check if the `"matrix"` node represents
//!      a valid R-vine matrix.
inline Vinecop::Vinecop(const boost::property_tree::ptree input,
                        const bool check_matrix)
{
    auto order =
        tools_serialization::ptree_to_vector<size_t>(input.get_child("order"));
    auto matrix =
        tools_serialization::ptree_to_triangular_array<size_t>(input.get_child("matrix"));

    vine_struct_ = RVineStructure(order, matrix, check_matrix);
    d_ = static_cast<size_t>(vine_struct_.get_dim());

    boost::property_tree::ptree pcs_node = input.get_child("pair copulas");
    for (size_t tree = 0; tree < d_ - 1; ++tree) {
        boost::property_tree::ptree tree_node;
        try {
            tree_node = pcs_node.get_child("tree" + std::to_string(tree));
        } catch (...) {
            break;  // vine was truncated, no more trees to parse
        }
        // reserve space for pair copulas of this tree
        pair_copulas_.resize(tree + 1);
        pair_copulas_[tree].resize(d_ - tree - 1);

        for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
            boost::property_tree::ptree pc_node =
                tree_node.get_child("pc" + std::to_string(edge));
            pair_copulas_[tree][edge] = Bicop(pc_node);
        }
    }

    threshold_ = 0;
    loglik_ = NAN;
}

//! creates from a JSON file
//! @param filename the name of the JSON file to read (see to_ptree() for the
//! structure of the file).
//! @param check_matrix whether to check if the `"matrix"` node represents
//!      a valid R-vine matrix.
inline Vinecop::Vinecop(const char *filename, const bool check_matrix) :
    Vinecop(tools_serialization::json_to_ptree(filename), check_matrix)
{
}

//! constructs a vine copula model from data by creating a model and calling
//! select_family().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param vine_struct an RVineStructure object specifying the vine structure.
//! @param controls see FitControlsVinecop.
inline Vinecop::Vinecop(const Eigen::MatrixXd &data,
        const RVineStructure &vine_struct,
        FitControlsVinecop controls)
{
    d_ = data.cols();
    nobs_ = data.rows();
    if (data.rows() == 1) {
        throw std::runtime_error("data must have more than one row");
    }
    if (d_ != vine_struct.get_dim()) {
        throw std::runtime_error("data and structure have "
                                     "incompatible dimensions.");
    }
    vine_struct_ = vine_struct;
    select_families(data, controls);
}

//! constructs a vine copula model from data by creating a model and calling
//! select_family().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param matrix either an empty matrix (default) or an R-vine structure
//!     matrix, see select_families().
//! @param controls see FitControlsVinecop.
//! @param check_matrix whether to check if `matrix` is a valid R-vine
//!     matrix.
inline Vinecop::Vinecop(const Eigen::MatrixXd &data,
        const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
        FitControlsVinecop controls,
        const bool check_matrix) :
    Vinecop(data,
            RVineStructure(matrix, check_matrix),
            controls) {}

//! constructs a vine copula model from data by creating a model and calling
//! select_family().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param order the order of the variables in the vine structure, see
//! the corresponding constructor of RVineStructure.
//! @param struct_array a triangular array object specifying the vine structure,
//! see the corresponding constructor of RVineStructure.
//! @param check_array whether `order` and `struct_array` shall be checked
//! for validity.
inline Vinecop::Vinecop(const Eigen::MatrixXd &data,
                        const std::vector<size_t> &order,
                        const TriangularArray<size_t> &struct_array,
                        FitControlsVinecop controls,
                        const bool check_array) :
    Vinecop(data,
            RVineStructure(order, struct_array, false, check_array),
            controls) {}

//! constructs a vine copula model from data by creating a model and
//! calling select_all().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param controls see FitControlsVinecop.
inline Vinecop::Vinecop(const Eigen::MatrixXd &data,
                        const FitControlsVinecop &controls)
{
    d_ = data.cols();
    nobs_ = data.rows();
    if (data.rows() == 1) {
        throw std::runtime_error("data must have more than one row");
    }
    select_all(data, controls);
}

//! Convert the copula into a boost::property_tree::ptree object, which
//! contains two nodes : `"matrix"` and `"pair copulas"`.
//! The former is encodes the R-Vine structure and the latter is a list of
//! child nodes for the trees (`"tree1"`, `"tree2"`, etc), each containing
//! a list of child nodes for the edges (`"pc1"`, `"pc2"`, etc).
//! See Bicop::to_ptree() for the encoding of pair-copulas.
//!
//! @return the boost::property_tree::ptree object containing the copula.
inline boost::property_tree::ptree Vinecop::to_ptree() const
{
    boost::property_tree::ptree pair_copulas;
    for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
        boost::property_tree::ptree tree_node;
        for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
            tree_node.add_child("pc" + std::to_string(edge),
                                pair_copulas_[tree][edge].to_ptree());
        }
        pair_copulas.add_child("tree" + std::to_string(tree), tree_node);
    }

    boost::property_tree::ptree output;
    output.add_child("pair copulas", pair_copulas);
    auto matrix_node = tools_serialization::triangular_array_to_ptree(get_struct_array());
    output.add_child("matrix", matrix_node);
    auto order_node = tools_serialization::vector_to_ptree(get_order());
    output.add_child("order", order_node);

    return output;
}

//! Write the copula object into a JSON file. See to_ptree() for the structure
//! of the file.
//!
//! @param filename the name of the file to write.
inline void Vinecop::to_json(const char *filename) const
{
    boost::property_tree::write_json(filename, to_ptree());
}

//! Initialize object for storing pair copulas
//!
//! @param d dimension of the vine copula.
//! @param truncation_level a truncation level (optional).
//! @return A nested vector such that `pc_store[t][e]` contains a Bicop.
//!     object for the pair copula corresponding to tree `t` and edge `e`.
inline std::vector<std::vector<Bicop>> Vinecop::make_pair_copula_store(
    const size_t d,
    const size_t truncation_level)
{
    return tools_select::VinecopSelector::make_pair_copula_store(
        d, truncation_level);
}


//! automatically fits and selects a vine copula model. Selection of the
//! structure is performed using the algorithm of
//! Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
//! *Selecting and estimating regular vine copulae and application to
//! financial returns.* Computational Statistics & Data Analysis, 59 (1),
//! 52-69.
//!
//! @param data nxd matrix of copula data.
//! @param controls the controls to the algorithm (see FitControlsVinecop).
inline void Vinecop::select_all(const Eigen::MatrixXd &data,
                                const FitControlsVinecop &controls)
{
    tools_eigen::check_if_in_unit_cube(data);
    check_data_dim(data);

    tools_select::StructureSelector selector(data, controls);
    if (controls.needs_sparse_select()) {
        selector.sparse_select_all_trees(data);
    } else {
        selector.select_all_trees(data);
    }
    threshold_ = selector.get_threshold();
    loglik_ = selector.get_loglik();
    nobs_ = data.rows();
    vine_struct_ = selector.get_rvine_matrix();
    pair_copulas_ = selector.get_pair_copulas();
}

//! automatically selects all pair-copula families and fits all parameters.
//!
//! @param data nxd matrix of copula data.
//! @param controls the controls to the algorithm (see FitControlsVinecop).
inline void Vinecop::select_families(const Eigen::MatrixXd &data,
                                     const FitControlsVinecop &controls)
{
    tools_eigen::check_if_in_unit_cube(data);
    check_data_dim(data);

    if (vine_struct_.get_trunc_lvl() > 0) {
        auto revorder = vine_struct_.get_order();
        tools_stl::reverse(revorder);
        auto newdata = data;
        for (size_t j = 0; j < d_; ++j)
            newdata.col(j) = data.col(revorder[j] - 1);

        tools_select::FamilySelector selector(newdata, vine_struct_, controls);
        if (controls.needs_sparse_select()) {
            selector.sparse_select_all_trees(newdata);
        } else {
            selector.select_all_trees(newdata);
        }
        threshold_ = selector.get_threshold();
        loglik_ = selector.get_loglik();
        nobs_ = data.rows();
        pair_copulas_ = selector.get_pair_copulas();
    }
}

//! @name Getters
//! @{

//! extracts a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline Bicop Vinecop::get_pair_copula(const size_t tree, const size_t edge) const
{
    if (tree > d_ - 2) {
        std::stringstream message;
        message <<
                "tree index out of bounds" << std::endl <<
                "allowed: 0, ..., " << d_ - 2 << std::endl <<
                "actual: " << tree << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    if (edge > d_ - tree - 2) {
        std::stringstream message;
        message <<
                "edge index out of bounds" << std::endl <<
                "allowed: 0, ..., " << d_ - tree - 2 << std::endl <<
                "actual: " << edge << std::endl <<
                "tree level: " << tree << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    if (tree >= pair_copulas_.size()) {
        // vine is truncated
        return Bicop();
    }
    return pair_copulas_[tree][edge];
}

//! extracts all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<Bicop>> Vinecop::get_all_pair_copulas() const
{
    return pair_copulas_;
}

//! extracts the family of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline BicopFamily Vinecop::get_family(const size_t tree, const size_t edge) const
{
    return get_pair_copula(tree, edge).get_family();
}

//! extracts the families of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<BicopFamily>> Vinecop::get_all_families() const
{
    std::vector<std::vector<BicopFamily>> families(pair_copulas_.size());
    for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
        families[tree].resize(d_ - 1 - tree);
        for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
            families[tree][edge] = get_family(tree, edge);
        }
    }

    return families;
}

//! extracts the rotation of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline int Vinecop::get_rotation(const size_t tree, const size_t edge) const
{
    return get_pair_copula(tree, edge).get_rotation();
}

//! extracts the rotations of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<int>> Vinecop::get_all_rotations() const
{
    std::vector<std::vector<int>> rotations(pair_copulas_.size());
    for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
        rotations[tree].resize(d_ - 1 - tree);
        for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
            rotations[tree][edge] = get_rotation(tree, edge);
        }
    }

    return rotations;
}

//! extracts the parameters of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline Eigen::MatrixXd Vinecop::get_parameters(const size_t tree, const size_t edge) const
{
    return get_pair_copula(tree, edge).get_parameters();
}

//! extracts the Kendall's \f$ tau \f$ of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline double Vinecop::get_tau(const size_t tree, const size_t edge) const
{
    return get_pair_copula(tree, edge).get_tau();
}


//! extracts the parameters of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<Eigen::MatrixXd>>
Vinecop::get_all_parameters() const
{
    std::vector<std::vector<Eigen::MatrixXd>>
        parameters(pair_copulas_.size());
    for (size_t tree = 0; tree < parameters.size(); ++tree) {
        parameters[tree].resize(d_ - 1 - tree);
        for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
            parameters[tree][edge] = get_parameters(tree, edge);
        }
    }

    return parameters;
}

//! extracts the Kendall's \f$ tau \f$s of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<double>> Vinecop::get_all_taus() const
{
    std::vector<std::vector<double>> taus(pair_copulas_.size());
    for (size_t tree = 0; tree < taus.size(); ++tree) {
        taus[tree].resize(d_ - 1 - tree);
        for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
            taus[tree][edge] = get_tau(tree, edge);
        }
    }

    return taus;
}

//! extracts the order vector of the vine copula model.
inline std::vector<size_t>
Vinecop::get_order() const
{
    return vine_struct_.get_order();
}

//! extracts the structure matrix of the vine copula model.
inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
Vinecop::get_matrix() const
{
    return vine_struct_.get_matrix();
}

//! extracts the above diagonal coefficients of the vine copula model.
inline TriangularArray<size_t>
Vinecop::get_struct_array() const
{
    return vine_struct_.get_struct_array();
}

//! extracts the log-likelihood (throws an error if model has not been 
//! fitted to data).
inline double Vinecop::get_loglik() const
{
    if (std::isnan(loglik_)) {
        throw std::runtime_error("copula has not been fitted from data ");
    }
    return loglik_;
}

//! extracts the number of observations used for the fit (throws an error if 
//! model has not been fitted to data).
inline size_t Vinecop::get_nobs() const
{
    if (std::isnan(loglik_)) {
        throw std::runtime_error("copula has not been fitted from data ");
    }
    return nobs_;
}

//! extracts the AIC (throws an error if model has not been 
//! fitted to data).
inline double Vinecop::get_aic() const
{
    if (std::isnan(loglik_)) {
        throw std::runtime_error("copula has not been fitted from data ");
    }
    return -2 * loglik_ + 2 * calculate_npars();
}

//! extracts the BIC (throws an error if model has not been 
//! fitted to data).
inline double Vinecop::get_bic() const
{
    if (std::isnan(loglik_)) {
        throw std::runtime_error("copula has not been fitted from data ");
    }
    return -2 * loglik_ + calculate_npars() * log(nobs_);;
}

//! extracts the log-likelihood (throws an error if model has not been 
//! fitted to data).
inline double Vinecop::get_mbicv(const double pi) const
{
    if (std::isnan(loglik_)) {
        throw std::runtime_error("copula has not been fitted from data ");
    }
    return -2 * loglik_ + this->compute_mbicv_penalty(nobs_, pi);
}

//! computes the penalty term for mBICV
inline double Vinecop::compute_mbicv_penalty(const size_t nobs, 
                                             const double pi) const
{
    if (!(pi > 0.0) | !(pi < 1.0)) {
        throw std::runtime_error("pi must be in the interval (0, 1)");
    }    
    auto all_fams = get_all_families();
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> non_indeps(d_ - 1);
    non_indeps.setZero();
    for (size_t t = 0; t < d_ - 1; t++) {
        if (t == all_fams.size()) {
            break;
        }
        for (size_t e = 0; e < d_ - 1 - t; e++) {
            if (all_fams[t][e] != BicopFamily::indep) {
                non_indeps(t)++;
            }
        }
    }
    auto sq0 = tools_stl::seq_int(1, d_ - 1);
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> sq(d_ - 1);
    auto pis = Eigen::VectorXd(d_ - 1);
    for (size_t i = 0; i < d_ - 1; i++) {
        sq(i) = sq0[i];
        pis(i) = std::pow(pi, sq0[i]);
    }
    double npars = this->calculate_npars();
    double log_prior = (
        non_indeps.cast<double>().array() * pis.array().log() +
        (d_ - non_indeps.array() - sq.array()).cast<double>() * 
        (1 - pis.array()).log()
    ).sum();
    
    return std::log(nobs) * npars - 2 * log_prior;
}

//! extracts the threshold (usually zero except `select_threshold == TRUE` in
//! `FitControlsVinecop()`).
inline double Vinecop::get_threshold() const
{
    return threshold_;
}


//! @}

//! calculates the density function of the vine copula model.
//!
//! @param u \f$ n \times d \f$ matrix of evaluation points.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::VectorXd Vinecop::pdf(const Eigen::MatrixXd &u, 
                                    const size_t num_threads) const
{
    tools_eigen::check_if_in_unit_cube(u);
    check_data_dim(u);
    size_t d = u.cols();
    size_t n = u.rows();

    // info about the vine structure (reverse rows (!) for more natural indexing)
    size_t trunc_lvl = pair_copulas_.size();
    std::vector<size_t> revorder;
    TriangularArray<size_t> no_matrix, max_matrix, needed_hfunc1, needed_hfunc2;
    if (trunc_lvl > 0) {
        revorder = vine_struct_.get_order();
        tools_stl::reverse(revorder);
        no_matrix = vine_struct_.get_struct_array();
        max_matrix = vine_struct_.get_max_array();
        needed_hfunc1 = vine_struct_.get_needed_hfunc1();
        needed_hfunc2 = vine_struct_.get_needed_hfunc2();
    }

    // initial value must be 1.0 for multiplication
    Eigen::VectorXd vine_density = Eigen::VectorXd::Constant(u.rows(), 1.0);

    auto do_batch = [&](const tools_batch::Batch& b) {

        // temporary storage objects for h-functions
        Eigen::MatrixXd hfunc1(b.size, d);
        Eigen::MatrixXd hfunc2(b.size, d);
        Eigen::MatrixXd u_e(b.size, 2);

        // fill first row of hfunc2 matrix with evaluation points;
        // points have to be reordered to correspond to natural order
        for (size_t j = 0; j < d; ++j)
            hfunc2.col(j) = u.block(b.begin, revorder[j] - 1, b.size, 1);

        for (size_t tree = 0; tree < trunc_lvl; ++tree) {
            tools_interface::check_user_interrupt(n * d > 1e5);
            for (size_t edge = 0; edge < d - tree - 1; ++edge) {
                tools_interface::check_user_interrupt(edge % 100 == 0);
                // extract evaluation point from hfunction matrices (have been
                // computed in previous tree level)
                size_t m = max_matrix(tree, edge);
                u_e.col(0) = hfunc2.col(edge);
                if (m == no_matrix(tree, edge)) {
                    u_e.col(1) = hfunc2.col(d - m);
                } else {
                    u_e.col(1) = hfunc1.col(d - m);
                }

                Bicop edge_copula = get_pair_copula(tree, edge);
                vine_density.segment(b.begin, b.size) =
                    vine_density.segment(b.begin, b.size).cwiseProduct(edge_copula.pdf(u_e));

                // h-functions are only evaluated if needed in next step
                if (needed_hfunc1(tree, edge)) {
                    hfunc1.col(edge) = edge_copula.hfunc1(u_e);
                }
                if (needed_hfunc2(tree, edge)) {
                    hfunc2.col(edge) = edge_copula.hfunc2(u_e);
                }
            }
        }
    };

    if (trunc_lvl > 0) {
        tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
        pool.map(do_batch, tools_batch::create_batches(n, num_threads));
        pool.join();
    }

    return vine_density;
}

//! calculates the cumulative distribution of the vine copula model.
//!
//! @param u \f$ n \times d \f$ matrix of evaluation points.
//! @param N integer for the number of quasi-random numbers to draw
//! to evaluate the distribution (default: 1e4).
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will generate `n` samples concurrently in 
//!   `num_threads` batches.
inline Eigen::VectorXd
Vinecop::cdf(const Eigen::MatrixXd &u, const size_t N, const size_t num_threads) const
{
    if (d_ > 21201) {
        std::stringstream message;
        message << "cumulative distribution available for models of " <<
                "dimension 21201 or less. This model's dimension: " << d_
                << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    tools_eigen::check_if_in_unit_cube(u);
    check_data_dim(u);

    // Simulate N quasi-random numbers from the vine model
    Eigen::MatrixXd U(N, d_);
    if (d_ > 300) {
        U = tools_stats::sobol(N, d_);
    } else {
        U = tools_stats::ghalton(N, d_);
    }
    U = inverse_rosenblatt(U, num_threads);

    // Alternative: simulate N pseudo-random numbers from the vine model
    //auto U = simulate(N);

    size_t n = u.rows();
    Eigen::VectorXd vine_distribution(n);
    Eigen::ArrayXXd x(N, 1);
    Eigen::RowVectorXd temp(d_);
    for (size_t i = 0; i < n; i++) {
        tools_interface::check_user_interrupt(i % 1000 == 0);
        temp = u.block(i, 0, 1, d_);
        x = (U.rowwise() - temp).rowwise().maxCoeff().array();
        vine_distribution(i) = (x <= 0.0).count();
    }
    return vine_distribution / static_cast<double>(N);
}

//! simulates from a vine copula model, see inverse_rosenblatt().
//!
//! @param n number of observations.
//! @param qrng set to true for quasi-random numbers.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will generate `n` samples concurrently in 
//!   `num_threads` batches.
//! @param seeds seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @return An \f$ n \times d \f$ matrix of samples from the copula model.
inline Eigen::MatrixXd Vinecop::simulate(const size_t n, 
                                         const bool qrng,
                                         const size_t num_threads,
                                         const std::vector<int>& seeds) const
{
    Eigen::MatrixXd U(n, d_);
    if (qrng) {
        if (d_ > 300) {
            U = tools_stats::sobol(n, d_);
        } else {
            U = tools_stats::ghalton(n, d_);
        }
    } else {
        U = tools_stats::simulate_uniform(n, d_, seeds);
    }

    return inverse_rosenblatt(U, num_threads);
}

//! calculates the log-likelihood, which is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, ..., U_{d, i}), \f]
//! where \f$ c \f$ is the copula density pdf().
//!
//! @param u \f$n \times d\f$ matrix of observations.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double Vinecop::loglik(const Eigen::MatrixXd &u, 
                              const size_t num_threads) const
{
    return pdf(u, num_threads).array().log().sum();
}

//! calculates the Akaike information criterion (AIC), whic is defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! calculate_npars(). The AIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param u \f$n \times 2\f$ matrix of observations.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double Vinecop::aic(const Eigen::MatrixXd &u, 
                           const size_t num_threads) const
{
    return -2 * loglik(u, num_threads) + 2 * calculate_npars();
}

//! calculates the Bayesian information criterion (BIC), which is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! calculate_npars(). The BIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param u \f$n \times 2\f$ matrix of observations.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double Vinecop::bic(const Eigen::MatrixXd &u, 
                           const size_t num_threads) const
{
    return -2 * loglik(u, num_threads) + 
        calculate_npars() * log(static_cast<double>(u.rows()));
}

//! calculates the modified Bayesian information criterion for vines (mBICV), 
//! which is defined as
//! \f[ \mathrm{mBICV} = -2\, \mathrm{loglik} +  \ln(n) \nu, - 2 * 
//! \sum_{t=1}^(d - 1) \{q_t log(\psi_0^t) - (d - t - q_t) log(1 -\psi_0^t)\}\f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood, \f$ \nu \f$ is the
//! (effective) number of parameters of the model, \f$ t \f$ is the tree level 
//! \f$ \psi_0 \f$ is the prior probability of having a non-independence copula 
//! in the first tree, and \f$ q_t \f$ is the number of non-independence copulas
//! in tree \f$ t \f$; The vBIC is a consistent model selection criterion for 
//! parametric sparse vine copula models when \f$ d = o(\sqrt{n \ln n})\f$.
//!
//! @param u \f$n \times 2\f$ matrix of observations.
//! @param pi baseline prior probability of a non-independence copula.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double Vinecop::mbicv(const Eigen::MatrixXd &u, 
                             const double pi, 
                             const size_t num_threads) const
{
    
    double n = static_cast<double>(u.rows());
    double ll = this->loglik(u, num_threads);
    return -2 * ll + this->compute_mbicv_penalty(n, pi);;       
}

//! returns sum of the number of parameters for all pair copulas (see
//! Bicop::calculate_npars()).
inline double Vinecop::calculate_npars() const
{
    double npars = 0.0;
    for (auto &tree : pair_copulas_) {
        for (auto &pc : tree) {
            npars += pc.calculate_npars();
        }
    }
    return npars;
}


//! calculates the inverse Rosenblatt transform for a vine copula model,
//! which can be used for simulation: the
//! function applied to independent uniform variates resembles simulated
//! data from the vine copula model.
//!
//! If the problem is too large, it is split recursively into halves (w.r.t.
//! n, the number of observations).
//! "Too large" means that the required memory will exceed 1 GB. An
//! examplary configuration requiring less than 1 GB is \f$ n = 1000 \f$,
//! \f$d = 200\f$.
//!
//! @param u \f$ n \times d \f$ matrix of evaluation points.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::MatrixXd
Vinecop::inverse_rosenblatt(const Eigen::MatrixXd &u,
                            const size_t num_threads) const
{
    tools_eigen::check_if_in_unit_cube(u);
    check_data_dim(u);
    size_t n = u.rows();
    if (n < 1) {
        throw std::runtime_error("n must be at least one");
    }
    size_t d = u.cols();

    Eigen::MatrixXd U_vine = u;  // output matrix

    //                   (direct + indirect)    (U_vine)       (info matrices)
    size_t bytes_required = (8 * 2 * n * d * d) + (8 * n * d) + (4 * 4 * d * d);
    // if the problem is too large (requires more than 1 GB memory), split
    // the data into two halves and call simulate on the reduced data.
    if ((n > 1) & (bytes_required > 1e9)) {
        size_t n_half = n / 2;
        size_t n_left = n - n_half;
        U_vine.block(0, 0, n_half, d) =
            inverse_rosenblatt(u.block(0, 0, n_half, d));
        U_vine.block(n_half, 0, n_left, d) =
            inverse_rosenblatt(u.block(n_half, 0, n_left, d));
        return U_vine;
    }

    // info about the vine structure (in upper triangular matrix notation)
    size_t trunc_lvl = pair_copulas_.size();
    std::vector<size_t> revorder, inverse_order;
    TriangularArray<size_t> no_matrix, max_matrix, needed_hfunc1, needed_hfunc2;
    if (trunc_lvl > 0) {
        revorder = vine_struct_.get_order();
        tools_stl::reverse(revorder);
        inverse_order = tools_stl::invert_permutation(revorder);
        no_matrix = vine_struct_.get_struct_array();
        max_matrix = vine_struct_.get_max_array();
        needed_hfunc1 = vine_struct_.get_needed_hfunc1();
        needed_hfunc2 = vine_struct_.get_needed_hfunc2();
    }

    auto do_batch = [&](const tools_batch::Batch& b) {
        if (d > 2) {
            // temporary storage objects for (inverse) h-functions
            TriangularArray<Eigen::VectorXd> hinv2(d + 1, trunc_lvl + 1);
            TriangularArray<Eigen::VectorXd> hfunc1(d + 1, trunc_lvl + 1);
        
            // initialize with independent uniforms (corresponding to natural
            // order)
            for (size_t j = 0; j < d; ++j) {
                hinv2(std::min(trunc_lvl, d - j - 1), j) = 
                    u.block(b.begin, revorder[j] - 1, b.size, 1);
            }
            hfunc1(0, d - 1) = hinv2(0, d - 1);
        
            // loop through variables (0 is just the initial uniform)
            for (ptrdiff_t var = d - 2; var >= 0; --var) {
            
                tools_interface::check_user_interrupt(n * d > 1e5);            
                size_t tree_start = std::min(trunc_lvl - 1, d - var - 2);
                for (ptrdiff_t tree = tree_start; tree >= 0; --tree) {
                    Bicop edge_copula = get_pair_copula(tree, var);
                
                    // extract data for conditional pair
                    Eigen::MatrixXd U_e(b.size, 2);
                    size_t m = max_matrix(tree, var);
                    U_e.col(0) = hinv2(tree + 1, var);
                    if (m == no_matrix(tree, var)) {
                        U_e.col(1) = hinv2(tree, d - m);
                    } else {
                        U_e.col(1) = hfunc1(tree, d - m);
                    }
                
                    // inverse Rosenblatt transform simulates data for conditional pair
                    hinv2(tree, var) = edge_copula.hinv2(U_e);
                
                    // if required at later stage, also calculate hfunc2
                    if (var < static_cast<ptrdiff_t>(d_) - 1) {
                        if (needed_hfunc1(tree, var)) {
                            U_e.col(0) = hinv2(tree, var);
                            hfunc1(tree + 1, var) = edge_copula.hfunc1(U_e);
                        }
                    }
                }
            }
            // go back to original order
            for (size_t j = 0; j < d; j++) {
                U_vine.block(b.begin, j, b.size, 1) = hinv2(0, inverse_order[j]);
            }
        } else {
            U_vine.block(b.begin, 0, b.size, 1) = get_pair_copula(0, 0).hinv2(u);
        }
    };

    if (trunc_lvl > 0) {
        tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
        pool.map(do_batch, tools_batch::create_batches(n, num_threads));
        pool.join();
    }

    return U_vine;
}

//! checks if dimension d of the data matches the dimension of the vine.
inline void Vinecop::check_data_dim(const Eigen::MatrixXd &data) const
{
    size_t d_data = data.cols();
    if (d_data != d_) {
        std::stringstream message;
        message << "data has wrong number of columns; " <<
                "expected: " << d_ <<
                ", actual: " << d_data << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
}

//! checks if pair copulas are compatible with the R-vine structure
inline  void Vinecop::check_pair_copulas_rvine_structure(
    const std::vector<std::vector<Bicop>> &pair_copulas) const
{
    size_t trunc_lvl = vine_struct_.get_trunc_lvl();
    if (pair_copulas.size() > std::min(d_ - 1, trunc_lvl)) {
        std::stringstream message;
        message << "pair_copulas is too large; "
                << "expected size: < " << std::min(d_ - 1, trunc_lvl) << ", "
                << "actual size: " << pair_copulas.size() << std::endl;
        throw std::runtime_error(message.str().c_str());
    }
    for (size_t t = 0; t < pair_copulas.size(); ++t) {
        if (pair_copulas[t].size() != d_ - 1 - t) {
            std::stringstream message;
            message << "size of pair_copulas[" << t << "] "
                    << "does not match dimension of matrix ("
                    << d_ << "); " << "expected size: "
                    << d_ - 1 - t << ", "
                    << "actual size: " << pair_copulas[t].size() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }
}

}
