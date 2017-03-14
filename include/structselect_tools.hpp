// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/graph_utility.hpp>
#include "vinecop_class.hpp"

// to allow for (auto e : boost::edges(g)) notation
namespace std
{
    template <class T>
    T begin(const std::pair<T,T>& eItPair) { return eItPair.first; }

    template <class T>
    T end(const std::pair<T,T>& eItPair) { return eItPair.second; }
}

namespace structselect_tools {
    // boost::graph represenation of a vine tree ----------------
    struct VertexProperties {
        std::vector<int> conditioning;
        std::vector<int> conditioned;
        std::vector<int> prev_edge_indices;
        VecXd hfunc1;
        VecXd hfunc2;
    };
    struct EdgeProperties {
        std::vector<int> conditioning;
        std::vector<int> conditioned;
        std::vector<int> all_indices;
        MatXd pc_data;
        VecXd hfunc1;
        VecXd hfunc2;
        double empirical_tau;
        BicopPtr pair_copula;
    };
    typedef boost::adjacency_list <
        boost::vecS,
        boost::vecS,
        boost::undirectedS,
        VertexProperties,
        boost::property<boost::edge_weight_t, double, EdgeProperties>
    > VineTree;


    // functions for manipulation of trees ----------------
    VineTree make_base_tree(const MatXd& data);
    VineTree select_next_tree(
        VineTree& prev_tree,
        std::vector<int> family_set,
        std::string selection_criterion,
        std::string method,
        bool preselect_families
    );
    VineTree edges_as_vertices(const VineTree& prev_tree);
    void add_allowed_edges(VineTree& tree);
    int find_common_neighbor(int v0, int v1, const VineTree& tree);
    MatXd get_pc_data(int v0, int v1, const VineTree& tree);
    void min_spanning_tree(VineTree &tree);
    void add_edge_info(VineTree& tree);
    void remove_edge_data(VineTree& tree);
    void remove_vertex_data(VineTree& tree);
    void select_pair_copulas(
        VineTree& tree,
        std::vector<int> family_set,
        std::string method,
        std::string selection_criterion,
        bool preselect_families
    );
    Vinecop as_vinecop(std::vector<VineTree>& trees);
    void flip(BicopPtr& bicop);
    void print_pair_copulas(VineTree& tree);
    std::string get_pc_index(
        boost::graph_traits<VineTree>::edge_descriptor e,
        VineTree& tree
    );


    // inline utility functions ----------------

    inline double pairwise_ktau(MatXd& u)
    {
        double tau;
        int n = u.rows();
        int two = 2;
        ktau_matrix(u.data(), &two, &n, &tau);
        return tau;
    }
}
