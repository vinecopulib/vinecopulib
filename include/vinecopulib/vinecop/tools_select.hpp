// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <cstddef>
#include <boost/graph/adjacency_list.hpp>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/vinecop/rvine_matrix.hpp>

// to allow for (auto e : boost::edges(g)) notation
namespace std
{
    template <class T>
    T begin(const std::pair<T,T>& eItPair) { return eItPair.first; }

    template <class T>
    T end(const std::pair<T,T>& eItPair) { return eItPair.second; }
}
namespace vinecopulib {

namespace tools_select {

double calculate_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data, 
                           std::string tree_criterion);
Eigen::MatrixXd calculate_criterion_matrix(const Eigen::MatrixXd& data, 
                                           std::string tree_criterion);
double calculate_gic(double loglik, double npars, int n);

class VinecopSelector 
{
public:
    VinecopSelector() {}
    ~VinecopSelector() {}
    
    virtual void select_next_tree() = 0;
    virtual void show_trace() = 0;
    void truncate();

protected:
    size_t n_;
    size_t d_;
    size_t trees_fitted_;
    FitControlsVinecop controls_;
};

class FamilySelector : public VinecopSelector 
{
public:
    FamilySelector(const Eigen::MatrixXd& data, 
                   const RVineMatrix& vine_matrix,
                   const FitControlsVinecop& controls);
    ~FamilySelector() {}
    
    void select_next_tree();
    void show_trace() {}  // TODO
    std::vector<std::vector<Bicop>> get_pair_copulas() {return pair_copulas_;}
    
private:
    // info about the vine structure
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> order_;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> no_matrix_;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> max_matrix_;
    tools_eigen::MatrixXb needed_hfunc1_;
    tools_eigen::MatrixXb needed_hfunc2_;
    // temporary storage objects for h-functions
    Eigen::MatrixXd hfunc1_;
    Eigen::MatrixXd hfunc2_;
    std::vector<std::vector<Bicop>> pair_copulas_;  // fitted pair copulas
};


// boost::graph represenation of a vine tree ----------------
struct VertexProperties {
    std::vector<size_t> conditioning;
    std::vector<size_t> conditioned;
    std::vector<size_t> prev_edge_indices;
    Eigen::VectorXd hfunc1;
    Eigen::VectorXd hfunc2;
};
struct EdgeProperties {
    std::vector<size_t> conditioning;
    std::vector<size_t> conditioned;
    std::vector<size_t> all_indices;
    Eigen::Matrix<double, Eigen::Dynamic, 2> pc_data;
    Eigen::VectorXd hfunc1;
    Eigen::VectorXd hfunc2;
    double weight;
    double crit;
    vinecopulib::Bicop pair_copula;
    double loglik;
    double npars;
    double fit_id;
};
typedef boost::adjacency_list <
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    VertexProperties,
    boost::property<boost::edge_weight_t, double, EdgeProperties>
> VineTree;

typedef std::pair<
boost::graph_traits<VineTree>::edge_descriptor, 
bool
> FoundEdge;

class StructureSelector : public VinecopSelector
{
public:
    StructureSelector(const Eigen::MatrixXd& data, 
                      const FitControlsVinecop& controls);
    ~StructureSelector() {}
    void select_next_tree();
    void show_trace();
    std::vector<VineTree> get_vine_trees() {return trees_;}
    
private:
    // functions for manipulation of trees ----------------
    VineTree make_base_tree(const Eigen::MatrixXd& data);
    VineTree select_next_tree(VineTree& prev_tree,
                              vinecopulib::FitControlsVinecop& controls,
                              const VineTree& tree_opt = VineTree());
    VineTree edges_as_vertices(const VineTree& prev_tree);
    void add_allowed_edges(VineTree& tree, std::string tree_criterion,
                           double threshold);
    ptrdiff_t find_common_neighbor(size_t v0, size_t v1, const VineTree& tree);
    Eigen::MatrixXd get_pc_data(size_t v0, size_t v1, const VineTree& tree);
    void min_spanning_tree(VineTree &tree);
    void add_edge_info(VineTree& tree);
    void remove_edge_data(VineTree& tree);
    void remove_vertex_data(VineTree& tree);
    void select_pair_copulas(VineTree& tree,
                             vinecopulib::FitControlsVinecop& controls);

    void select_pair_copulas(VineTree& tree,
        vinecopulib::FitControlsVinecop& controls,
        const VineTree& tree_opt);
    FoundEdge find_old_fit(double fit_id, const VineTree& old_graph);
    std::vector<double> get_thresholded_edge_crits(
        const std::vector<VineTree>& trees,
        vinecopulib::FitControlsVinecop& controls);
    double get_next_threshold(std::vector<double>& thresholded_crits,
                              double learning_rate = 0.025);
    double get_tree_loglik(const VineTree& tree);
    double get_tree_npars(const VineTree& tree);

    void print_pair_copulas(VineTree& tree);
    std::string get_pc_index(boost::graph_traits<VineTree>::edge_descriptor e,
                             VineTree& tree);
    
    std::vector<VineTree> trees_;
};

}

}
