// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

// to allow for (auto e : boost::edges(g)) notation
namespace std {
template<class T>
T
begin(const std::pair<T, T>& eItPair)
{
  return eItPair.first;
}

template<class T>
T
end(const std::pair<T, T>& eItPair)
{
  return eItPair.second;
}
}
namespace vinecopulib {

namespace tools_select {

double
calculate_criterion(const Eigen::MatrixXd& data,
                    std::string tree_criterion,
                    Eigen::VectorXd weights);

Eigen::MatrixXd
calculate_criterion_matrix(const Eigen::MatrixXd& data,
                           const std::string& tree_criterion,
                           const Eigen::VectorXd& weights);

std::vector<size_t>
get_disc_cols(std::vector<std::string> var_types);

// boost::graph represenation of a vine tree
struct VertexProperties
{
  std::vector<size_t> conditioning;
  std::vector<size_t> conditioned;
  std::vector<size_t> all_indices;
  std::vector<size_t> prev_edge_indices;
  Eigen::VectorXd hfunc1;
  Eigen::VectorXd hfunc2;
  Eigen::VectorXd hfunc1_sub;
  Eigen::VectorXd hfunc2_sub;
  std::vector<std::string> var_types{ "c", "c" };
};
struct EdgeProperties
{
  std::vector<size_t> conditioning;
  std::vector<size_t> conditioned;
  std::vector<size_t> all_indices;
  Eigen::MatrixXd pc_data;
  Eigen::VectorXd hfunc1;
  Eigen::VectorXd hfunc2;
  Eigen::VectorXd hfunc1_sub;
  Eigen::VectorXd hfunc2_sub;
  std::vector<std::string> var_types{ "c", "c" };
  double weight;
  double crit;
  vinecopulib::Bicop pair_copula;
  double fit_id;
};
typedef boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  VertexProperties,
  boost::property<boost::edge_weight_t, double, EdgeProperties>>
  VineTree;

typedef boost::graph_traits<VineTree>::edge_descriptor EdgeIterator;
typedef std::pair<EdgeIterator, bool> FoundEdge;
typedef boost::property_map<VineTree, boost::edge_weight_t>::type WeightMap;

class VinecopSelector
{
public:
  VinecopSelector(const Eigen::MatrixXd& data,
                  const FitControlsVinecop& controls,
                  std::vector<std::string> var_types);

  VinecopSelector(const Eigen::MatrixXd& data,
                  const RVineStructure& vine_struct,
                  const FitControlsVinecop& controls,
                  std::vector<std::string> var_types);

  virtual ~VinecopSelector() = default;

  std::vector<std::vector<Bicop>> get_pair_copulas() const;

  RVineStructure get_rvine_structure() const;

  static std::vector<std::vector<Bicop>> make_pair_copula_store(
    size_t d,
    size_t trunc_lvl);

  void select_all_trees(const Eigen::MatrixXd& data);

  void sparse_select_all_trees(const Eigen::MatrixXd& data);

  double get_loglik() const;

  double get_threshold() const;

  size_t get_nobs() const;

  std::vector<VineTree> get_trees() const { return trees_; };
  std::vector<VineTree> get_trees_opt() const { return trees_opt_; };

protected:
  virtual void select_tree(size_t t);

  void finalize(size_t trunc_lvl);

  double get_mbicv_of_tree(size_t t, double loglik);

  double get_loglik_of_tree(size_t t);

  double get_npars_of_tree(size_t t);

  void set_tree_to_indep(size_t t);

  void print_pair_copulas_of_tree(size_t);

  std::vector<double> get_thresholded_crits();

  void initialize_new_fit(const Eigen::MatrixXd& data);

  void set_current_fit_as_opt(const double& loglik);

  void add_allowed_edges(VineTree& vine_tree);

  void select_edges(VineTree& vine_tree);

  Eigen::MatrixXd get_pc_data(size_t v0, size_t v1, const VineTree& tree);

  Eigen::VectorXd get_hfunc(const VertexProperties& vertex_data, bool is_first);

  Eigen::VectorXd get_hfunc_sub(const VertexProperties& vertex_data,
                                bool is_first);

  ptrdiff_t find_common_neighbor(size_t v0, size_t v1, const VineTree& tree);

  virtual double compute_fit_id(const EdgeProperties& e);

  size_t n_;
  size_t d_;
  bool structure_unknown_{ true };
  std::vector<std::string> var_types_;
  FitControlsVinecop controls_;
  tools_thread::ThreadPool pool_;
  std::vector<VineTree> trees_;
  RVineStructure vine_struct_;
  std::vector<std::vector<Bicop>> pair_copulas_;
  // for sparse selction
  std::vector<VineTree> trees_opt_;
  double loglik_;
  double threshold_;
  double psi0_; // initial prior probability for mbicv

  double get_next_threshold(std::vector<double>& thresholded_crits);

  // functions for manipulation of trees ----------------
  VineTree make_base_tree(const Eigen::MatrixXd& data);

  VineTree edges_as_vertices(const VineTree& prev_tree);

  void min_spanning_tree(VineTree& tree);

  void add_edge_info(VineTree& tree);

  void add_pc_info(const EdgeIterator& e, VineTree& tree);

  void remove_edge_data(VineTree& tree);

  void remove_vertex_data(VineTree& tree);

  void select_pair_copulas(VineTree& tree,
                           const VineTree& tree_opt = VineTree());

  FoundEdge find_old_fit(double fit_id, const VineTree& old_graph);

  double get_tree_loglik(const VineTree& tree);

  double get_tree_npars(const VineTree& tree);

  size_t get_num_non_indeps_of_tree(size_t t);

  std::string get_pc_index(const EdgeIterator& e, const VineTree& tree);
};
}
}

#include <vinecopulib/vinecop/implementation/tools_select.ipp>
