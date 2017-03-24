// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "vinecop/tools_structselect.hpp"
#include "misc/tools_stl.hpp"
#include "misc/tools_stats.hpp"
#include <iostream>
#include <cmath>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace tools_structselect {
    
    using namespace tools_stl;

    //! Create base tree of the vine
    //!
    //!  The base tree is a star on d + 1 variables, where the conditioning
    //!  set of each edge consists of a single number. When building the next
    //!  tree, the edges become vertices. Because the base graph was a star
    //!  all edges are allowed by the proximity condition, and the edges will
    //!  have a conditioning set consisting of the two vertex indices. This
    //!  will be the first actual tree of the vine.
    //!
    //!  @param data nxd matrix of copula data.
    //!  @return A VineTree object containing the base graph.
    VineTree make_base_tree(const Eigen::MatrixXd& data)
    {
        int d = data.cols();
        VineTree base_tree(d);
        // a star connects the root node (d) with all other nodes
        for (int target = 0; target < d; ++target) {
            // add edge and extract edge iterator
            auto e = add_edge(d, target, base_tree).first;

            // inititialize hfunc1 with actual data for variable "target"
            base_tree[e].hfunc1 = data.col(boost::target(e, base_tree));
            // identify edge with variable "target" and initialize sets
            base_tree[e].conditioning.reserve(2);
            base_tree[e].conditioning.push_back(boost::target(e, base_tree));
            base_tree[e].conditioned.reserve(d - 2);
        }

        return base_tree;
    }

    //! Select and fit next tree of the vine
    //!
    //! The next tree is found the following way:
    //!     1. Edges of the previous tree become edges in the new tree.
    //!     2. All edges allowed by the proximity condition are added to the new
    //!        graph.
    //!     3. Collapse the new graph to a maximum spanning tree for edge 
    //!        weight.
    //!     4. Populate edges with conditioning/conditioned sets and pseudo-
    //!        observations.
    //!     5. Fit and select a copula model for each edge.
    //!
    //! @param prev_tree tree T_{k}.
    //! @param family_set the set of copula families to consider (if empty, then
    //!     all families are included; all families are included by default).
    //! @param method indicates the estimation method: either maximum likelihood
    //!     estimation (method = "mle", default) or inversion of Kendall's tau
    //!     (method = "itau"). When method = "itau" is used with families having
    //!     more thanone parameter, the main dependence parameter is found by
    //!     inverting the Kendall's tau and the remainders by profile likelihood
    //!     optimization.
    //! @param threshold for thresholded vines.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param selection_criterion the selection criterion; either "aic" or "bic"
    //!     (default).
    //! @param preselect_families  whether to exclude families before fitting based
    //!     on symmetry properties of the data.
    //! @return tree T_{k+1}.
    VineTree select_next_tree(
        VineTree& prev_tree,
        std::vector<vinecopulib::BicopFamily> family_set,
        std::string method,
        double threshold,
        std::string tree_criterion,
        std::string selection_criterion,
        bool preselect_families
    )
    {
        auto new_tree = edges_as_vertices(prev_tree);
        remove_edge_data(prev_tree); // no longer needed
        add_allowed_edges(new_tree, tree_criterion, threshold);
        if (boost::num_vertices(new_tree) > 2) {
            min_spanning_tree(new_tree);
        }
        add_edge_info(new_tree);  // for pc estimation and next tree
        remove_vertex_data(new_tree);  // no longer needed
        select_pair_copulas(
            new_tree,
            family_set,
            method,
            threshold,
            selection_criterion,
            preselect_families
        );

        return new_tree;
    }

    //! Convert edge set into vertex set of a new graph
    //!
    //! Further information about the structure is passed along:
    //!     - conditioning/conditioned set,
    //!     - indices of vertices connected by the edge in the previous tree.
    //!
    //! @param tree T_{k}.
    //! @return A edge-less graph of vertices, each representing one edge of the
    //!     previous tree.
    VineTree edges_as_vertices(const VineTree& prev_tree)
    {
        // start with full graph
        int d = num_edges(prev_tree);
        VineTree new_tree(d);

        // cut & paste information from previous tree
        int i = 0;
        for (auto e : boost::edges(prev_tree)) {
            new_tree[i].hfunc1 = prev_tree[e].hfunc1;
            new_tree[i].hfunc2 = prev_tree[e].hfunc2;
            new_tree[i].conditioning = prev_tree[e].conditioning;
            new_tree[i].conditioned = prev_tree[e].conditioned;
            new_tree[i].prev_edge_indices.reserve(2);
            new_tree[i].prev_edge_indices.push_back(boost::source(e, prev_tree));
            new_tree[i].prev_edge_indices.push_back(boost::target(e, prev_tree));
            ++i;
        }

        return new_tree;
    }

    //! Add edges allowed by the proximity condition
    //!
    //! Also calculates the edge weight (e.g., 1-|tau| for tree_criterion = 
    //! "itau").
    //!
    //! @param vine_tree tree of a vine.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param threshold for thresholded vines.
    void add_allowed_edges(
        VineTree& vine_tree, 
        std::string tree_criterion,
        double threshold
    )
    {
        for (auto v0 : boost::vertices(vine_tree)) {
            for (unsigned int v1 = 0; v1 < v0; ++v1) {
                // check proximity condition: common neighbor in previous tree
                // (-1 means 'no common neighbor')
                if (find_common_neighbor(v0, v1, vine_tree) > -1) {
                    auto pc_data = get_pc_data(v0, v1, vine_tree);
                    auto w = get_tree_criterion(pc_data, tree_criterion, threshold);
                    auto e = boost::add_edge(v0, v1, w, vine_tree).first;
                    vine_tree[e].weight = w;
                }
            }
        }
    }
        
    double get_tree_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data, 
        std::string tree_criterion, double threshold) 
    {
        double w;
        if (tree_criterion == "tau") {
            w = 1.0 - std::fabs(tools_stats::pairwise_ktau(data));
        } else if (tree_criterion == "hoeffd") {
            // scale to [0,1]
            w = 1.0 - (30*tools_stats::pairwise_hoeffd(data)+0.5)/1.5;
        } else if (tree_criterion == "rho") {
            w = 1.0 - std::fabs(tools_stats::pairwise_cor(data));
        } else {
            throw std::runtime_error("tree criterion not implemented");
        }
        if (w > 1 - threshold) {
            w = 1.0;
        }
        
        return w;
    }

    // Find common neighbor in previous tree
    //
    // @param v0,v1 vertices in the tree.
    // @param tree the current tree.
    // @return Gives the index of the vertex in the previous tree that was
    //     shared by e0, e1, the edge representations of v0, v1.
    int find_common_neighbor(int v0, int v1, const VineTree& tree)
    {
        auto ei0 = tree[v0].prev_edge_indices;
        auto ei1 = tree[v1].prev_edge_indices;
        auto ei_common = intersect(ei0, ei1);

        if (ei_common.size() == 0) {
            return -1;
        } else {
            return ei_common[0];
        }
    }

    // Extract pair copula pseudo-observations from h-functions
    //
    // @param v0,v1 vertex indices.
    // @param tree a vine tree.
    // @return The pseudo-observations for the pair coula, extracted from
    //     the h-functions calculated in the previous tree.
    Eigen::MatrixXd get_pc_data(int v0, int v1, const VineTree& tree)
    {
        Eigen::MatrixXd pc_data(tree[v0].hfunc1.size(), 2);
        int ei_common = find_common_neighbor(v0, v1, tree);
        if (find_position(ei_common, tree[v0].prev_edge_indices) == 0) {
            pc_data.col(0) = tree[v0].hfunc1;
        } else {
            pc_data.col(0) = tree[v0].hfunc2;
        }
        if (find_position(ei_common, tree[v1].prev_edge_indices) == 0) {
            pc_data.col(1) = tree[v1].hfunc1;
        } else {
            pc_data.col(1) = tree[v1].hfunc2;
        }

        return pc_data;
    }

    //! Collapse a graph to the minimum spanning tree
    //!
    //! @param graph the input graph.
    //! @return the input graph with all non-MST edges removed.
    void min_spanning_tree(VineTree &graph)
    {
        int d =  num_vertices(graph);
        std::vector<int> targets(d);
        prim_minimum_spanning_tree(graph, targets.data());
        for (int v1 = 0; v1 < d; ++v1) {
            for (int v2 = 0; v2 < v1; ++v2) {
                if ((v2 != targets[v1]) & (v1 != targets[v2])) {
                    boost::remove_edge(v1, v2, graph);
                }
            }
        }
    }

    //! Add conditioning info and data for each edge
    //!
    //! See, e.g., Czado (2010), "Pair-copula constructions of multivariate
    //! copulas", url: https://mediatum.ub.tum.de/doc/1079253/file.pdf
    //! @param tree a vine tree.
    void add_edge_info(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            auto v0 = boost::source(e, tree);
            auto v1 = boost::target(e, tree);
            tree[e].pc_data = get_pc_data(v0, v1, tree);

            auto v0_indices = cat(tree[v0].conditioning, tree[v0].conditioned);
            auto v1_indices = cat(tree[v1].conditioning, tree[v1].conditioned);

            auto test = intersect(v0_indices, v1_indices);
            auto d01 = set_diff(v0_indices, v1_indices);
            auto d10 = set_diff(v1_indices, v0_indices);

            tree[e].conditioning = cat(d01, d10);
            tree[e].conditioned = intersect(v0_indices, v1_indices);
            tree[e].all_indices = cat(tree[e].conditioning, tree[e].conditioned);
        }
    }

    //! Remove data (hfunc1/hfunc2/pc_data) from all edges of a vine tree
    //! @param tree a vine tree.
    void remove_edge_data(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            tree[e].hfunc1 = Eigen::VectorXd();
            tree[e].hfunc2 = Eigen::VectorXd();
            tree[e].pc_data = Eigen::Matrix<double, Eigen::Dynamic, 2>(0, 2);
        }
    }

    //! Remove data (hfunc1/hfunc2/pc_data) from all vertices of a vine tree
    //! @param tree a vine tree.
    void remove_vertex_data(VineTree& tree)
    {
        for (auto v : boost::vertices(tree)) {
            tree[v].hfunc1 = Eigen::VectorXd();
            tree[v].hfunc2 = Eigen::VectorXd();
        }
    }

    //! Fit and select a pair copula for each edges
    //! @param tree a vine tree preprocessed with add_edge_info().
    //! @param family_set the set of copula families to consider (if empty, then
    //!     all families are included; all families are included by default).
    //! @param method indi::cates the estimation method: either maximum likelihood
    //!     estimation (method = "mle", default) or inversion of Kendall's tau
    //!     (method = "itau"). When method = "itau" is used with families having
    //!     more thanone parameter, the main dependence parameter is found by
    //!     inverting the Kendall's tau and the remainders by profile likelihood
    //!     optimization.
    //! @param the threshold for thresholded vines.
    //! @param selection_criterion the selection criterion; either "aic" or "bic"
    //!     (default).
    //! @param preselect_families  whether to exclude families before fitting based
    //!     on symmetry properties of the data.
    void select_pair_copulas(
        VineTree& tree,
        std::vector<vinecopulib::BicopFamily> family_set,
        std::string method,
        double threshold,
        std::string selection_criterion,
        bool preselect_families
    )
    {
        for (auto e : boost::edges(tree)) {
            if (tree[e].weight > 1 - threshold) {
                tree[e].pair_copula = vinecopulib::Bicop();
            } else {
                tree[e].pair_copula = vinecopulib::Bicop(
                    tree[e].pc_data,
                    family_set,
                    method,
                    selection_criterion,
                    preselect_families
                );
            }

            tree[e].hfunc1 = tree[e].pair_copula.hfunc1(tree[e].pc_data);
            tree[e].hfunc2 = tree[e].pair_copula.hfunc2(tree[e].pc_data);
        }
    }

    //! Convert fitted trees into Vinecop object
    //!
    //! @param trees a vector of trees preprocessed by add_edge_info(); the
    //!     0th entry should be the base graph and is not used.
    //! @return Vinecop object corresponding to the fitted trees.
    vinecopulib::Vinecop as_vinecop(std::vector<VineTree>& trees)
    {
        int d = trees.size();
        Eigen::MatrixXi mat = Eigen::MatrixXi::Constant(d, d, 0);
        auto pcs = vinecopulib::Vinecop::make_pair_copula_store(d);

        for (int col = 0; col < d - 1; ++col) {
            int t = d - 1 - col;
            // start with highest tree in this column and fill first two
            // entries by conditioning set
            auto e0 = *boost::edges(trees[t]).first;
            mat(t, col) = trees[t][e0].conditioning[0];
            mat(t - 1, col) = trees[t][e0].conditioning[1];

            // assign fitted pair copula to appropriate entry, see
            // vinecopulib::Vinecop::get_pair_copula().
            pcs[t - 1][col] = trees[t][e0].pair_copula;

            // initialize running set with full conditioing set of this edge
            auto ned_set = trees[t][e0].conditioned;

            // iteratively search for an edge in lower tree that shares all indices
            // in the conditioning set + diagonal entry
            for (int k = 1; k < t; ++k) {
                auto reduced_set = cat(mat(t, col), ned_set);
                for (auto e : boost::edges(trees[t - k])) {
                    if (is_same_set(trees[t - k][e].all_indices, reduced_set)) {
                        // next matrix entry is conditioning variable of new edge
                        // that's not equal to the diagonal entry of this column
                        auto e_new = trees[t - k][e];
                        auto pos = find_position(mat(t, col), e_new.conditioning);
                        mat(t - k - 1, col) = e_new.conditioning[std::abs(1 - pos)];
                        if (pos == 1) {
                            e_new.pair_copula.flip();
                        }
                        // assign fitted pair copula to appropriate entry, see
                        // vinecopulib::Vinecop::get_pair_copula().
                        pcs[t - 1 - k][col] = e_new.pair_copula;

                        // start over with conditioning set of next edge
                        ned_set = e_new.conditioned;

                        // remove edge (must not be reused in another column!)
                        int v0 = boost::source(e, trees[t - k]);
                        int v1 = boost::target(e, trees[t - k]);
                        boost::remove_edge(v0, v1, trees[t - k]);
                        break;
                    }
                }
            }
        }

        // The last column contains a single element which must be different
        // from all other diagonal elements. Based on the properties of an
        // R-vine matrix, this must be the element next to it.
        mat(0, d - 1) = mat(0, d - 2);

        // change to user-facing format
        // (variable index starting at 1 instead of 0)
        Eigen::MatrixXi new_mat = mat;
        for (int i = 0; i < d; ++i)
            for (int j = 0; j < d - i; ++j)
                new_mat(i, j) += 1;

        return vinecopulib::Vinecop(pcs, new_mat);
    }
    
    //! Print indices, family, and parameters for each pair-copula
    //! @param tree a vine tree.
    void print_pair_copulas(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            std::stringstream pc_info;
            pc_info <<
                get_pc_index(e, tree) << " <-> " <<
                "fam = " << tree[e].pair_copula.get_family_name() <<
                ", rot = " << tree[e].pair_copula.get_rotation() <<
                ", par = " <<  tree[e].pair_copula.get_parameters() <<
                std::endl;
            std::cout << pc_info.str().c_str();
        }
    }

    //! Get edge index for the vine (like 1, 2; 3)
    //! @param e a descriptor for the edge.
    //! @param tree a vine tree.
    std::string get_pc_index(
        boost::graph_traits<VineTree>::edge_descriptor e,
        VineTree& tree
     )
     {
        std::stringstream index;
        // add 1 everywhere for user-facing representation (boost::graph starts
        // at 0)
        index <<
            tree[e].conditioning[0] + 1 <<
            "," <<
            tree[e].conditioning[1] + 1;
        if (tree[e].conditioned.size() > 0) {
            index << " ; ";
            for (unsigned int i = 0; i < tree[e].conditioned.size(); ++i) {
                index << tree[e].conditioned[i] + 1;
                if (i < tree[e].conditioned.size() - 1)
                    index << ",";
                }
        }

        return index.str().c_str();
    }

}
