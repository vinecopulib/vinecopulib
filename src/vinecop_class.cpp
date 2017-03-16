// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "vinecop_class.hpp"
#include "tools_structselect.hpp"
#include "bicop_indep.hpp"
#include "tools_stl.hpp"
#include "tools_stats.hpp"
#include <vector>
#include <exception>
#include <iostream>

namespace vinecopulib
{
    //! Construct a vine copula object of dimension d
    //!
    //! @param d dimension of the vine copula.
    //!
    //! @return a d-dimensional D-vine with variable order 1, ..., d and all
    //!     pair-copulas set to independence
    Vinecop::Vinecop(int d)
    {
        d_ = d;
        // D-vine with variable order (1, ..., d)
        Eigen::VectorXi order(d);
        for (int i = 0; i < d; ++i) {
            order(i) = i + 1;
        }
        vine_matrix_ = RVineMatrix(RVineMatrix::construct_d_vine_matrix(order));

        // all pair-copulas are independence
        pair_copulas_ = make_pair_copula_store(d);
        for (auto& tree : pair_copulas_) {
            for (auto& pc : tree) {
                pc = BicopPtr(new IndepBicop);
            }
        }
    }

    //! Construct a vine copula object from a vector<BicopPtr> and structure matrix
    //!
    //! @param pair_copulas a nested vector of BicopPtrs; can be initialized by
    //!     make_pair_copula_store(d).
    //! @param matrix R-vine matrix.
    //!
    //! @return a d-dimensional D-vine with variable order 1, ..., d and all
    //!     pair-copulas set to independence
    Vinecop::Vinecop(
            const std::vector<std::vector<BicopPtr>>& pair_copulas,
            const Eigen::MatrixXi& matrix
    )
    {
        d_ = matrix.rows();
        if ((int) pair_copulas.size() != d_ - 1) {
            std::stringstream message;
            message <<
                    "size of of pair_copulas does not match dimension of matrix (" <<
                    d_ << ")" <<
                    "expected size:" << d_ - 1 << ", "<<
                    "actual size:" << pair_copulas.size() << std::endl;
            throw std::runtime_error(message.str().c_str());

        }
        for (int t = 0; t < d_ - 1; ++t) {
            if ((int) pair_copulas[t].size() != d_ - 1 - t) {
                std::stringstream message;
                message <<
                        "size of of pair_copulas[" << t << "] " <<
                        "does not match dimension of matrix (" << d_ << ")" <<
                        "expected size:" << d_ - 1 - t << ", "<<
                        "actual size:" << pair_copulas[t].size() << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
        }

        // D-vine with variable order (1, ..., d)
        vine_matrix_ = RVineMatrix(matrix);
        pair_copulas_ = pair_copulas;
    }

    //! Initialize object for storing pair copulas
    //!
    //! @param d dimension of the vine copula.
    //! @return A nested vector such that pc_store(d)[t][e] contains a BicopPtr to
    //!     the pair copula corresponding to tree t and edge e.
    std::vector<std::vector<BicopPtr>> Vinecop::make_pair_copula_store(int d)
    {
        std::vector<std::vector<BicopPtr>> pc_store(d - 1);
        for (int t = 0; t < d - 1; ++t) {
            pc_store[t].resize(d - 1 - t);
        }

        return pc_store;
    }


    //! Automated model and structure selection for vine copulas
    //!
    //! Implements the structure selection algorithm of  Dissmann et al. (2013).
    //!
    //! @param data nxd matrix of copula data.
    //! @param family_set the set of copula families to consider (if empty, then
    //!     all families are included; all families are included by default).
    //! @param method indicates the estimation method: either maximum likelihood
    //!     estimation (method = "mle", default) or inversion of Kendall's tau
    //!     (method = "itau"). When method = "itau" is used with families having
    //!     more thanone parameter, the main dependence parameter is found by
    //!     inverting the Kendall's tau and the remainders by profile likelihood
    //!     optimization.
    //! @param selection_criterion the selection criterion; either "aic" or "bic"
    //!     (default).
    //! @param preselect_families  whether to exclude families before fitting based
    //!     on symmetry properties of the data.
    //! @param show_trace whether to show a trace of the building progress (default
    //!     is false).
    //! @return The fitted vine copula model.
    Vinecop Vinecop::select(
            const Eigen::MatrixXd& data,
            std::vector<BicopFamily> family_set,
            std::string method,
            int truncation_level,
            Eigen::MatrixXi matrix,
            std::string selection_criterion,
            bool preselect_families,
            bool show_trace
    )
    {
        using namespace tools_structselect;
        int d = data.cols();
        if (matrix.size() > 0) {
            throw std::runtime_error("fixed matrix selection not implemented yet.");
        }
        std::vector<VineTree> trees(d);

        trees[0] = make_base_tree(data);
        for (int t = 1; t < d; ++t) {
            // select tree structure and pair copulas
            trees[t] = select_next_tree(
                    trees[t - 1],
                    family_set,
                    method,
                    selection_criterion,
                    preselect_families
            );

            // print out fitted pair-copulas for this tree
            if (show_trace) {
                std::cout << "Tree " << t - 1 << ":" << std::endl;
                print_pair_copulas(trees[t]);
            }

            // truncate (only allow for Independence copula from here on)
            if (truncation_level == t) {
                family_set = std::vector<BicopFamily>({BicopFamily::Indep});
            }
        }

        return as_vinecop(trees);;
    }

    //! Access a pair copula
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    //!
    //! @return A \code std::shared_ptr (alias \code BicopPtr) to a \code Bicop
    //!     object.
    BicopPtr Vinecop::get_pair_copula(int tree, int edge)
    {
        if (tree > d_ - 2) {
            std::stringstream message;
            message <<
                    "tree index out of bounds" << std::endl <<
                    "allowed: 0, ..., " << d_ - 2 << std::endl <<
                    "actual: " << tree << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        if ((edge < 0) | (edge > d_ - tree - 2)) {
            std::stringstream message;
            message <<
                    "edge index out of bounds" << std::endl <<
                    "allowed: 0, ..., " << d_ - tree - 2 << std::endl <<
                    "actual: " << edge << std::endl <<
                    "tree level: " <<  tree  << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        return pair_copulas_[tree][edge];
    }

    //! Get family of a pair copula
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    //!
    //! @return An \code int containing the family index.
    BicopFamily Vinecop::get_family(int tree, int edge)
    {
        return get_pair_copula(tree, edge)->get_family();
    }
    
    //! Get families of all pair copulas
    //!
    //! @return A nested vector containing in \c vec[t][e] the family for
    //!     edge e in tree t.
    std::vector<std::vector<BicopFamily>> Vinecop::get_all_families()
    {
        std::vector<std::vector<BicopFamily>> families(d_ - 1);
        for (int t = 0; t < d_ - 1; ++t)
            families[t].resize(d_ - 1 - t);
        for (int tree = 0; tree < d_ - 1; ++tree) {
            for (int edge = 0; edge < d_ - 1 - tree; ++edge) {
                families[tree][edge] = get_family(tree, edge);
            }
        }

        return families;
    }

    //! Get rotation of a pair copula
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    //!
    //! @return An \code int containing the rotation.
    int Vinecop::get_rotation(int tree, int edge)
    {
        return get_pair_copula(tree, edge)->get_rotation();
    }

    //! Get rotations of all pair copulas
    //!
    //! @return A nested vector containing in \c vec[t][e] the rotation for
    //!     edge e in tree t.
    std::vector<std::vector<int>> Vinecop::get_all_rotations()
    {
        std::vector<std::vector<int>> rotations(d_ - 1);
        for (int t = 0; t < d_ - 1; ++t)
            rotations[t].resize(d_ - 1 - t);
        for (int tree = 0; tree < d_ - 1; ++tree) {
            for (int edge = 0; edge < d_ - 1 - tree; ++edge) {
                rotations[tree][edge] = get_rotation(tree, edge);
            }
        }

        return rotations;
    }

    //! Get parameters of a pair copula
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    //!
    //! @return An \code Eigen::VectorXd containing the parameters.
    Eigen::VectorXd Vinecop::get_parameters(int tree, int edge)
    {
        return get_pair_copula(tree, edge)->get_parameters();
    }

    //! Get parameters of all pair copulas
    //!
    //! @return A nested vector containing in \c vec[t][e] the parameters for
    //!     edge e in tree t.
    std::vector<std::vector<Eigen::VectorXd>> Vinecop::get_all_parameters()
    {
        std::vector<std::vector<Eigen::VectorXd>> parameters(d_ - 1);
        for (int t = 0; t < d_ - 1; ++t)
            parameters[t].resize(d_ - 1 - t);
        for (int tree = 0; tree < d_ - 1; ++tree) {
            for (int edge = 0; edge < d_ - 1 - tree; ++edge) {
                parameters[tree][edge] = get_parameters(tree, edge);
            }
        }

        return parameters;
    }
    
    //! Probability density function of a vine copula
    //!
    //! @param u mxd matrix of evaluation points.
    Eigen::VectorXd Vinecop::pdf(const Eigen::MatrixXd& u)
    {
        int d = u.cols();
        int n = u.rows();
        if (d != d_) {
            std::stringstream message;
            message << "u has wrong number of columns. " <<
                    "expected: " << d_ <<
                    ", actual: " << d << std::endl;
            throw std::runtime_error(message.str().c_str());
        }

        // info about the vine structure (reverse rows (!) for more natural indexing)
        Eigen::VectorXi revorder      = vine_matrix_.get_order().reverse();
        Eigen::MatrixXi no_matrix     = vine_matrix_.in_natural_order();
        Eigen::MatrixXi max_matrix    = vine_matrix_.get_max_matrix();
        MatrixXb needed_hfunc1 = vine_matrix_.get_needed_hfunc1();
        MatrixXb needed_hfunc2 = vine_matrix_.get_needed_hfunc2();

        // initial value must be 1.0 for multiplication
        Eigen::VectorXd vine_density = Eigen::VectorXd::Constant(u.rows(), 1.0);

        // temporary storage objects for h-functions
        Eigen::MatrixXd hfunc1(n, d);
        Eigen::MatrixXd hfunc2(n, d);
        Eigen::MatrixXd u_e(n, 2);

        // fill first row of hfunc2 matrix with evaluation points;
        // points have to be reordered to correspond to natural order
        for (int j = 0; j < d; ++j)
            hfunc2.col(j) = u.col(revorder(j) - 1);

        for (int tree = 0; tree < d - 1; ++tree) {
            for (int edge = 0; edge < d - tree - 1; ++edge) {
                // get pair copula for this edge
                BicopPtr edge_copula = get_pair_copula(tree, edge);

                // extract evaluation point from hfunction matrices (have been
                // computed in previous tree level)
                int m = max_matrix(tree, edge);
                u_e.col(0) = hfunc2.col(edge);
                if (m == no_matrix(tree, edge)) {
                    u_e.col(1) = hfunc2.col(d - m);
                } else {
                    u_e.col(1) = hfunc1.col(d - m);
                }

                vine_density = vine_density.cwiseProduct(edge_copula->pdf(u_e));
                // h-functions are only evaluated if needed in next step
                if (needed_hfunc1(tree + 1, edge)) {
                    hfunc1.col(edge) = edge_copula->hfunc1(u_e);
                }
                if (needed_hfunc2(tree + 1, edge)) {
                    hfunc2.col(edge) = edge_copula->hfunc2(u_e);
                }
            }
        }

        return vine_density;
    }

    //! Simulate from a vine copula model
    //!
    //! If the problem is too large, it is split recursively into halves (w.r.t
    //! to n, the number of observations).
    //! "Too large" means that the required memory will exceed 1 GB. An examplary
    //! configuration requiring less than 1GB is n = 1000, d = 200.
    //!
    //! @param n number of observations.
    //! @return Simulated data from the vine copula model.
    Eigen::MatrixXd Vinecop::simulate(int n)
    {
        Eigen::MatrixXd U = tools_stats::simulate_uniform(n, d_);
        return inverse_rosenblatt(U);
    }


    //! Inverse rosenblatt transform for a VineCopula model
    //! 
    //! If the problem is too large, it is split recursively into halves (w.r.t
    //! to n, the number of observations).
    //! "Too large" means that the required memory will exceed 1 GB. An examplary
    //! configuration requiring less than 1GB is n = 1000, d = 200.
    //! 
    //! @param U mxd matrix of indpendent uniform random variables.
    //! @return The transformed data corresponds to simulated data from the
    //!     vine copula model when U are independent unfiorms.
    Eigen::MatrixXd Vinecop::inverse_rosenblatt(const Eigen::MatrixXd& U)
    {
        int n = U.rows();
        if (n < 1) {
            throw std::runtime_error("n must be at least one");
        }
        int d = U.cols();
        if (d != d_) {
            std::stringstream message;
            message << "U has wrong number of columns; " <<
                    "expected: " << d_ <<
                    ", actual: " << d << std::endl;
            throw std::runtime_error(message.str().c_str());
        }

        Eigen::MatrixXd U_vine = U;  // output matrix

        //                   (direct + indirect)    (U_vine)       (info matrices)
        int bytes_required = (8 * 2 * n * d * d) +  (8 * n * d)  + (4 * 4 * d * d);
        // if the problem is too large (requires more than 1 GB memory), split
        // the data into two halves and call simulate on the reduced data.
        if ((n > 1) & (bytes_required > 1e9)) {
            int n_half = n / 2;
            int n_left = n - n_half;
            U_vine.block(0, 0, n_half, d) = 
                inverse_rosenblatt(U.block(0, 0, n_half, d));
            U_vine.block(n_half, 0, n_left, d) = 
                inverse_rosenblatt(U.block(n_half, 0, n_left, d));
            return U_vine;
        }

        // info about the vine structure (in upper triangular matrix notation)
        Eigen::VectorXi revorder      = vine_matrix_.get_order().reverse();
        Eigen::MatrixXi no_matrix     = vine_matrix_.in_natural_order();
        Eigen::MatrixXi max_matrix    = vine_matrix_.get_max_matrix();
        MatrixXb needed_hfunc1 = vine_matrix_.get_needed_hfunc1();
        MatrixXb needed_hfunc2 = vine_matrix_.get_needed_hfunc2();

        // temporary storage objects for (inverse) h-functions
        typedef Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> Array3d;
        Array3d hinv2(d, d);
        Array3d hfunc1(d, d);

        // initialize with independent uniforms (corresponding to natural order)
        for (int j = 0; j < d; ++j)
            hinv2(d - j - 1, j) = U.col(revorder(j) - 1);
        hfunc1(0, d - 1) = hinv2(0, d - 1);

        // loop through variables (0 is just the inital uniform)
        for (int var = d - 2; var > -1; --var) {
            for (int tree = d - var - 2; tree > -1; --tree) {
                BicopPtr edge_copula = get_pair_copula(tree, var);

                // extract data for conditional pair
                Eigen::MatrixXd U_e(n, 2);
                int m = max_matrix(tree, var);
                U_e.col(0) = hinv2(tree + 1, var);
                if (m == no_matrix(tree, var)) {
                    U_e.col(1) = hinv2(tree, d - m);
                } else {
                    U_e.col(1) = hfunc1(tree, d - m);
                }

                // inverse Rosenblatt transform simulates data for conditional pair
                hinv2(tree, var) = edge_copula->hinv2(U_e);

                // if required at later stage, also calculate hfunc2
                if (var < d_ - 1) {
                    if (needed_hfunc1(tree + 1, var)) {
                        U_e.col(0) = hinv2(tree, var);
                        hfunc1(tree + 1, var) = edge_copula->hfunc1(U_e);
                    }
                }
            }
        }

        // go back to original order
        Eigen::VectorXi inverse_order = inverse_permutation(revorder);
        for (int j = 0; j < d; ++j)
            U_vine.col(j) = hinv2(0, inverse_order(j));

        return U_vine;
    }

    // get indexes for reverting back to old order in simulation routine
    Eigen::VectorXi inverse_permutation(const Eigen::VectorXi& order) {
        // start with (0, 1, .., k)
        auto indexes = tools_stl::seq_int(0, order.size());

        // get sort indexes by comparing values in order
        std::sort(indexes.begin(), indexes.end(),
                  [&order](int i1, int i2) {return order(i1) < order(i2);});

        // convert to Eigen::VectorXi;
        return Eigen::Map<Eigen::VectorXi>(&indexes[0], order.size());
    }
}
