// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

namespace vinecopulib {
//! @brief A class for vine copula models
//!
//! A vine copula model is characterized by the structure matrix (see
//! RVineMatrix) and the pair-copulas.
class Vinecop
{
public:
    // Constructors
    Vinecop()
    {
    }

    Vinecop(size_t d);

    Vinecop(const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
            const bool check_matrix = true);

    Vinecop(const std::vector <std::vector<Bicop>> &pair_copulas,
            const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
            const bool check_matrix = true);

    Vinecop(const Eigen::MatrixXd &data,
            const FitControlsVinecop &controls = FitControlsVinecop());

    Vinecop(const Eigen::MatrixXd &data,
            const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
            FitControlsVinecop controls = FitControlsVinecop(),
            const bool check_matrix = true);

    Vinecop(const char *filename, const bool check_matrix = true);

    Vinecop(const boost::property_tree::ptree input,
            const bool check_matrix = true);

    // Serialize
    boost::property_tree::ptree to_ptree() const;

    void to_json(const char *filename) const;

    // Methods modifying structure and/or families and parameters
    void select_all(const Eigen::MatrixXd &data,
                    const FitControlsVinecop &controls = FitControlsVinecop());

    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop());

    // Getters for a single pair copula
    Bicop get_pair_copula(size_t tree, size_t edge) const;

    BicopFamily get_family(size_t tree, size_t edge) const;

    int get_rotation(size_t tree, size_t edge) const;

    Eigen::MatrixXd get_parameters(size_t tree, size_t edge) const;
    
    double get_tau(size_t tree, size_t edge) const;

    // Getters for all pair copulas
    std::vector <std::vector<Bicop>> get_all_pair_copulas() const;

    std::vector <std::vector<BicopFamily>> get_all_families() const;

    std::vector <std::vector<int>> get_all_rotations() const;

    std::vector <std::vector<Eigen::MatrixXd>> get_all_parameters() const;
    
    std::vector <std::vector<double>> get_all_taus() const;

    // Getters for the structure
    std::vector<size_t> get_order() const;

    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    RVineMatrix<size_t> get_struct_matrix() const;
    
    // getter for the threshold
    double get_threshold() const;

    // getter for the loglik
    double get_loglik() const;

    // Stats methods
    Eigen::VectorXd pdf(const Eigen::MatrixXd &u,
                        const size_t num_threads = 1) const;

    Eigen::VectorXd cdf(const Eigen::MatrixXd &u, 
                        const size_t N = 1e4,
                        const size_t num_threads = 1) const;

    Eigen::MatrixXd simulate(const size_t n, 
                             const bool qrng = false, 
                             const size_t num_threads = 1,
                             const std::vector<int>& seeds =  {1, 2, 3, 4}) const;

    Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd &u,
                                       const size_t num_threads = 1) const;

    // Fit statistics
    double calculate_npars() const;

    double loglik(const Eigen::MatrixXd &u) const;

    double aic(const Eigen::MatrixXd &u) const;

    double bic(const Eigen::MatrixXd &u) const;
    
    double mbicv(const Eigen::MatrixXd &u, const double pi) const;

    // Misc methods
    static std::vector <std::vector<Bicop>>
    make_pair_copula_store(const size_t d,
                           const size_t truncation_level = std::numeric_limits<size_t>::max());

private:
    size_t d_;
    RVineStructure vine_struct_;
    std::vector <std::vector<Bicop>> pair_copulas_;
    double threshold_;
    double loglik_;

    void check_data_dim(const Eigen::MatrixXd &data) const;
};

}

#include <vinecopulib/vinecop/implementation/class.ipp>
