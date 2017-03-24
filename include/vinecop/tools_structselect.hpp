// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include "vinecop/class.hpp"
#include "bicop/family.hpp"

// to allow for (auto e : boost::edges(g)) notation
namespace std
{
    template <class T>
    T begin(const std::pair<T,T>& eItPair) { return eItPair.first; }

    template <class T>
    T end(const std::pair<T,T>& eItPair) { return eItPair.second; }
}

namespace tools_structselect {
    // boost::graph represenation of a vine tree ----------------
    struct VertexProperties {
        std::vector<int> conditioning;
        std::vector<int> conditioned;
        std::vector<int> prev_edge_indices;
        Eigen::VectorXd hfunc1;
        Eigen::VectorXd hfunc2;
    };
    struct EdgeProperties {
        std::vector<int> conditioning;
        std::vector<int> conditioned;
        std::vector<int> all_indices;
        Eigen::Matrix<double, Eigen::Dynamic, 2> pc_data;
        Eigen::VectorXd hfunc1;
        Eigen::VectorXd hfunc2;
        double weight;
        vinecopulib::Bicop pair_copula;
    };
    typedef boost::adjacency_list <
        boost::vecS,
        boost::vecS,
        boost::undirectedS,
        VertexProperties,
        boost::property<boost::edge_weight_t, double, EdgeProperties>
    > VineTree;


    // functions for manipulation of trees ----------------
    VineTree make_base_tree(const Eigen::MatrixXd& data);
    VineTree select_next_tree(
        VineTree& prev_tree,
        std::vector<vinecopulib::BicopFamily> family_set,
        std::string method,
        std::string tree_criterion,
        std::string selection_criterion,
        bool preselect_families
    );
    VineTree edges_as_vertices(const VineTree& prev_tree);
    void add_allowed_edges(VineTree& tree, std::string tree_criterion);
    double get_edge_weight(
        Eigen::Matrix<double, Eigen::Dynamic, 2> data,
        std::string tree_criterion
     );
    int find_common_neighbor(int v0, int v1, const VineTree& tree);
    Eigen::MatrixXd get_pc_data(int v0, int v1, const VineTree& tree);
    void min_spanning_tree(VineTree &tree);
    void add_edge_info(VineTree& tree);
    void remove_edge_data(VineTree& tree);
    void remove_vertex_data(VineTree& tree);
    void select_pair_copulas(
        VineTree& tree,
        std::vector<vinecopulib::BicopFamily> family_set,
        std::string method,
        std::string selection_criterion,
        bool preselect_families
    );
    vinecopulib::Vinecop as_vinecop(std::vector<VineTree>& trees);
    //void flip(vinecopulib::Bicop& bicop);
    void print_pair_copulas(VineTree& tree);
    std::string get_pc_index(
        boost::graph_traits<VineTree>::edge_descriptor e,
        VineTree& tree
    );

}
