// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/vinecop/tools_select.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/vinecop/class.hpp>

#include <cmath>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace vinecopulib {
    
namespace tools_select {

using namespace tools_stl;

//! Calculate criterion for tree selection
//! @param data observations.
//! @param tree_criterion the criterion.
double calculate_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                           std::string tree_criterion)
{
    double w = 0.0;
    if (tree_criterion == "tau") {
        w = std::fabs(tools_stats::pairwise_ktau(data));
    } else if (tree_criterion == "hoeffd") {
        // scale to [0,1]
        w = (30 * tools_stats::pairwise_hoeffd(data) + 0.5) / 1.5;
    } else if (tree_criterion == "rho") {
        w = std::fabs(tools_stats::pairwise_cor(data));
    }
    
    return w;
}

//! Calculates maximal criterion for tree selection.
//! @param data observations.
//! @param tree_criterion the criterion.
Eigen::MatrixXd calculate_criterion_matrix(const Eigen::MatrixXd& data, 
                                           std::string tree_criterion)
{
    Eigen::MatrixXd w;
    if (tree_criterion == "tau") {
        w = tools_stats::ktau_matrix(data);
    } else if (tree_criterion == "hoeffd") {
        // scale to [0,1]
        w = (30 * tools_stats::hoeffd_matrix(data).array() + 0.5) / 1.5;
    } else if (tree_criterion == "rho") {
        w = tools_stats::cor_matrix(data);
    }

    return w.array().abs();
}

//! calculates the Generalized Information Criterion.
double calculate_gic(double loglik, double npars, int n)
{
    double log_npars = std::log(npars);
    if (npars == 0.0) {
        log_npars = 0.0;
    }
    return -2 * loglik + std::log(std::log(n)) * log_npars * npars;
}

std::vector<std::vector<Bicop>> VinecopSelector::get_pair_copulas() const
{
        return pair_copulas_;
}

RVineMatrix VinecopSelector::get_rvine_matrix() const
{
    return vine_matrix_;
}

//! truncates the vine copula at the current tree level (tracked by 
//! `trees_fitted_` data member).
void VinecopSelector::truncate() 
{
    controls_.set_family_set({BicopFamily::indep});
}

void VinecopSelector::select_all_trees()
{
    for (size_t t = 0; t < d_ - 1; ++t) {
        select_next_tree();  // select pair copulas (+ structure) of tree t
    
        if (controls_.get_show_trace()) {
            show_trace();  // print fitted pair-copulas for this tree
        }
        
        if (controls_.get_truncation_level() == t + 1) {
            truncate();  // only allow for Independence copula from here on
        }
    }
    finalize();
}

FamilySelector::FamilySelector(const Eigen::MatrixXd& data, 
                               const RVineMatrix& rvm,
                               const FitControlsVinecop& controls)
{
    n_ = data.rows();
    d_ = data.cols();
    order_ = rvm.get_order();
    no_matrix_ = rvm.in_natural_order();
    max_matrix_ = rvm.get_max_matrix();
    needed_hfunc1_ = rvm.get_needed_hfunc1();
    needed_hfunc2_ = rvm.get_needed_hfunc2();
    hfunc1_ = Eigen::MatrixXd::Zero(n_, d_);
    hfunc2_ = Eigen::MatrixXd::Zero(n_, d_);
    // fill first row of hfunc2 matrix with evaluation points;
    // points have to be reordered to correspond to natural order
    for (size_t j = 0; j < d_; ++j) {
        hfunc2_.col(j) = data.col(order_.reverse()(j) - 1);
    }
    trees_fitted_ = 0;
    controls_ = controls;
    vine_matrix_ = rvm;
    pair_copulas_ = Vinecop::make_pair_copula_store(d_);
}

void FamilySelector::select_next_tree()
{
    Eigen::MatrixXd u_e(n_, 2);
    size_t tree = trees_fitted_++;
    for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        size_t m = max_matrix_(tree, edge);
        u_e.col(0) = hfunc2_.col(edge);
        if (m == no_matrix_(tree, edge)) {
            u_e.col(1) = hfunc2_.col(d_ - m);
        } else {
            u_e.col(1) = hfunc1_.col(d_ - m);
        }
        
        // select pair-copula
        if (tree > controls_.get_truncation_level()) {
            pair_copulas_[tree][edge] = Bicop(BicopFamily::indep);
        } else {
            if (controls_.get_threshold() != 0) {
                double crit = 
                    calculate_criterion(u_e, controls_.get_tree_criterion());
                if (crit < controls_.get_threshold()) {
                    pair_copulas_[tree][edge] = Bicop(BicopFamily::indep);
                }
            } else {
                pair_copulas_[tree][edge] = Bicop(u_e, controls_);
            }
        }
        
        // h-functions are only evaluated if needed in next step
        if (needed_hfunc1_(tree + 1, edge)) {
            hfunc1_.col(edge) = pair_copulas_[tree][edge].hfunc1(u_e);
        }
        if (needed_hfunc2_(tree + 1, edge)) {
            hfunc2_.col(edge) = pair_copulas_[tree][edge].hfunc2(u_e);
        }
    }
}


StructureSelector::StructureSelector(const Eigen::MatrixXd& data, 
                                     const FitControlsVinecop& controls)
{
    n_ = data.rows();
    d_ = data.cols();
    trees_ = std::vector<VineTree>(d_);
    trees_[0] = make_base_tree(data);
    controls_ = controls;
    trees_fitted_ = 0;
    pair_copulas_ = Vinecop::make_pair_copula_store(d_);
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
//! @param controls the controls for fitting a vine copula 
//!     (see FitControlsVinecop).
//! @param tree_opt the current optimal tree (used only for sparse selection).
void StructureSelector::select_next_tree()
{
    auto new_tree = edges_as_vertices(trees_[trees_fitted_]);
    remove_edge_data(trees_[trees_fitted_]); // no longer needed
    add_allowed_edges(new_tree);
    if (boost::num_vertices(new_tree) > 2) {
        min_spanning_tree(new_tree);
    }
    add_edge_info(new_tree);       // for pc estimation and next tree
    remove_vertex_data(new_tree);  // no longer needed
    auto tree_opt = VineTree();    // TODO
    select_pair_copulas(new_tree, tree_opt);

    trees_[++trees_fitted_] =  new_tree;
}

//! prints pair copulas for recently fitted tree.
void StructureSelector::show_trace()
{
    std::stringstream tree_info;
    tree_info << "** Tree: " << trees_fitted_ << std::endl;
    tools_interface::print(tree_info.str().c_str());
    print_pair_copulas(trees_[trees_fitted_]);
}

void StructureSelector::finalize()
{
    using namespace tools_stl;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(d_, d_);
    mat.fill(0);

    for (size_t col = 0; col < d_ - 1; ++col) {
        size_t t = d_ - 1 - col;
        // start with highest tree in this column and fill first two
        // entries by conditioning set
        auto e0 = *boost::edges(trees_[t]).first;
        mat(t, col) = trees_[t][e0].conditioning[0];
        mat(t - 1, col) = trees_[t][e0].conditioning[1];

        // assign fitted pair copula to appropriate entry, see
        // vinecopulib::Vinecop::get_pair_copula().
        pair_copulas_[t - 1][col] = trees_[t][e0].pair_copula;

        // initialize running set with full conditioing set of this edge
        auto ned_set = trees_[t][e0].conditioned;

        // iteratively search for an edge in lower tree that shares all indices
        // in the conditioning set + diagonal entry
        for (size_t k = 1; k < t; ++k) {
            auto reduced_set = cat(mat(t, col), ned_set);
            for (auto e : boost::edges(trees_[t - k])) {
                if (is_same_set(trees_[t - k][e].all_indices, reduced_set)) {
                    // next matrix entry is conditioning variable of new edge
                    // that's not equal to the diagonal entry of this column
                    auto e_new = trees_[t - k][e];
                    ptrdiff_t pos = find_position(mat(t, col), e_new.conditioning);
                    mat(t - k - 1, col) = e_new.conditioning[std::abs(1 - pos)];
                    if (pos == 1) {
                        e_new.pair_copula.flip();
                    }
                    // assign fitted pair copula to appropriate entry, see
                    // vinecopulib::Vinecop::get_pair_copula().
                    pair_copulas_[t - 1 - k][col] = e_new.pair_copula;

                    // start over with conditioning set of next edge
                    ned_set = e_new.conditioned;

                    // remove edge (must not be reused in another column!)
                    size_t v0 = boost::source(e, trees_[t - k]);
                    size_t v1 = boost::target(e, trees_[t - k]);
                    boost::remove_edge(v0, v1, trees_[t - k]);
                    break;
                }
            }
        }
    }

    // The last column contains a single element which must be different
    // from all other diagonal elements. Based on the properties of an
    // R-vine matrix, this must be the element next to it.
    mat(0, d_ - 1) = mat(0, d_ - 2);

    // change to user-facing format
    // (variable index starting at 1 instead of 0)
    auto new_mat = mat;
    for (size_t i = 0; i < d_; ++i) {
        for (size_t j = 0; j < d_ - i; ++j) {
            new_mat(i, j) += 1;    
        }
    }
    vine_matrix_ = RVineMatrix(new_mat);
}

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
VineTree StructureSelector::make_base_tree(const Eigen::MatrixXd& data)
{
    size_t d = data.cols();
    VineTree base_tree(d);
    // a star connects the root node (d) with all other nodes
    for (size_t target = 0; target < d; ++target) {
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

//! Convert edge set into vertex set of a new graph
//!
//! Further information about the structure is passed along:
//!     - conditioning/conditioned set,
//!     - indices of vertices connected by the edge in the previous tree.
//!
//! @param tree T_{k}.
//! @return A edge-less graph of vertices, each representing one edge of the
//!     previous tree.
VineTree StructureSelector::edges_as_vertices(const VineTree& prev_tree)
{
    // start with full graph
    size_t d = num_edges(prev_tree);
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
void StructureSelector::add_allowed_edges(VineTree& vine_tree)
{
    std::string tree_criterion = controls_.get_tree_criterion();
    double threshold = controls_.get_threshold();
    for (auto v0 : boost::vertices(vine_tree)) {
        for (size_t v1 = 0; v1 < v0; ++v1) {
            // check proximity condition: common neighbor in previous tree
            // (-1 means 'no common neighbor')
            if (find_common_neighbor(v0, v1, vine_tree) > -1) {
                auto pc_data = get_pc_data(v0, v1, vine_tree);
                double crit = calculate_criterion(pc_data, tree_criterion);
                double w = 1.0 - (double)(crit > threshold) * crit;
                auto e = boost::add_edge(v0, v1, w, vine_tree).first;
                vine_tree[e].weight = w;
                vine_tree[e].crit = crit;
            }
        }
    }
}

// Find common neighbor in previous tree
//
// @param v0,v1 vertices in the tree.
// @param tree the current tree.
// @return Gives the index of the vertex in the previous tree that was
//     shared by e0, e1, the edge representations of v0, v1.
ptrdiff_t StructureSelector::find_common_neighbor(size_t v0, size_t v1, 
                                                  const VineTree& tree)
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
Eigen::MatrixXd StructureSelector::get_pc_data(size_t v0, size_t v1, 
                                               const VineTree& tree)
{
    Eigen::MatrixXd pc_data(tree[v0].hfunc1.size(), 2);
    size_t ei_common = find_common_neighbor(v0, v1, tree);
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
void StructureSelector::min_spanning_tree(VineTree &graph)
{
    size_t d =  num_vertices(graph);
    std::vector<size_t> targets(d);
    prim_minimum_spanning_tree(graph, targets.data());
    for (size_t v1 = 0; v1 < d; ++v1) {
        for (size_t v2 = 0; v2 < v1; ++v2) {
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
void StructureSelector::add_edge_info(VineTree& tree)
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
void StructureSelector::remove_edge_data(VineTree& tree)
{
    for (auto e : boost::edges(tree)) {
        tree[e].hfunc1 = Eigen::VectorXd();
        tree[e].hfunc2 = Eigen::VectorXd();
        tree[e].pc_data = Eigen::Matrix<double, Eigen::Dynamic, 2>(0, 2);
    }
}

//! Remove data (hfunc1/hfunc2/pc_data) from all vertices of a vine tree
//! @param tree a vine tree.
void StructureSelector::remove_vertex_data(VineTree& tree)
{
    for (auto v : boost::vertices(tree)) {
        tree[v].hfunc1 = Eigen::VectorXd();
        tree[v].hfunc2 = Eigen::VectorXd();
    }
}

//! Fit and select a pair copula for each edges
//! @param tree a vine tree preprocessed with add_edge_info().
//! @param tree_opt the current optimal tree (used only for sparse 
//!     selection).
//! @{
void StructureSelector::select_pair_copulas(VineTree& tree)
{
    for (auto e : boost::edges(tree)) {
        if (tree[e].crit < controls_.get_threshold()) {
            tree[e].pair_copula = vinecopulib::Bicop();
        } else {
            tree[e].pair_copula = vinecopulib::Bicop(tree[e].pc_data, controls_);
        }

        tree[e].hfunc1 = tree[e].pair_copula.hfunc1(tree[e].pc_data);
        tree[e].hfunc2 = tree[e].pair_copula.hfunc2(tree[e].pc_data);
    }
}

void StructureSelector::select_pair_copulas(VineTree& tree,
                                            const VineTree& tree_opt)
{
    for (auto e : boost::edges(tree)) {
        bool is_thresholded = (tree[e].crit < controls_.get_threshold());
        bool used_old_fit = false;
        // the formula is quite arbitrary, but sufficient for 
        // identifying situations where fits can be re-used
        tree[e].fit_id = tree[e].pc_data(0, 0) - 2 * tree[e].pc_data(0, 1); 
        tree[e].fit_id += 5.0 * (double) is_thresholded;
        if (boost::num_edges(tree_opt) > 0) {
            auto old_fit = find_old_fit(tree[e].fit_id, tree_opt);
            if (old_fit.second)  {  // second indicates if match was found
                // data and thresholding status haven't changed, 
                // we can use old fit
                used_old_fit = true;
                tree[e].pair_copula = tree_opt[old_fit.first].pair_copula;
            }
        }
        if (!used_old_fit) {
            if (is_thresholded) {
                tree[e].pair_copula = vinecopulib::Bicop();
            } else {
                tree[e].pair_copula.select(tree[e].pc_data, controls_);
            }
        }
        
        tree[e].hfunc1 = tree[e].pair_copula.hfunc1(tree[e].pc_data);
        tree[e].hfunc2 = tree[e].pair_copula.hfunc2(tree[e].pc_data);
        tree[e].loglik = tree[e].pair_copula.loglik(tree[e].pc_data);
        tree[e].npars  = tree[e].pair_copula.calculate_npars();
    }
}
//! @}


//! extracts all criterion values that got thresholded to zero.
std::vector<double> StructureSelector::thresholded_crits(
    const std::vector<VineTree>& trees)
{
    std::vector<double> out;
    for (size_t t = 1; t < trees.size(); ++t) {
        for (auto e : boost::edges(trees[t])) {
            if (trees[t][e].crit < controls_.get_threshold()) {
                out.push_back(trees[t][e].crit);
            }
        }
    }
    
    return out;
}

//! chooses threshold for next iteration such that at a proportion of at 
//! least `learning_rate` of the previously thresholded pairs become 
//! non-thresholded.
double StructureSelector::get_next_threshold(
    std::vector<double>& thresholded_crits,
    double learning_rate)
{
    if (thresholded_crits.size() == 0) {
        return 0.0;
    }
    // sort in descending order
    std::sort(thresholded_crits.begin(), thresholded_crits.end());
    std::reverse(thresholded_crits.begin(), thresholded_crits.end());
    size_t m = thresholded_crits.size();
    // pick threshold that changes at least <rate>% of the pair-copulas
    return thresholded_crits[std::ceil(m * learning_rate) - 1];
}

//! finds the fitted pair-copula from the previous iteration.
FoundEdge StructureSelector::find_old_fit(double fit_id, 
                                          const VineTree& old_graph) 
{
    auto edge = boost::edge(0, 1, old_graph).first;
    bool fit_with_same_id = false;
    for (auto e : boost::edges(old_graph)) {
        if (fit_id == old_graph[e].fit_id) {
            fit_with_same_id = true;
            edge = e;
        }
    }
    return std::make_pair(edge, fit_with_same_id);
}

//! calculates the log-likelihood of a tree.
double StructureSelector::get_tree_loglik(const VineTree& tree)
{
    double ll = 0.0;
    for (const auto& e : boost::edges(tree)) {
        ll += tree[e].loglik;
    }
    return ll;
}

//! calculates the numbers of parameters of a tree.
double StructureSelector::get_tree_npars(const VineTree& tree)
{
    double npars = 0.0;
    for (const auto& e : boost::edges(tree)) {
        npars += tree[e].npars;
    }
    return npars;
}

//! Print indices, family, and parameters for each pair-copula
//! @param tree a vine tree.
void StructureSelector::print_pair_copulas(VineTree& tree)
{
    for (auto e : boost::edges(tree)) {
        std::stringstream pc_info;
        pc_info << get_pc_index(e, tree) << " <-> " <<
                    tree[e].pair_copula.str() << std::endl;
        vinecopulib::tools_interface::print(pc_info.str().c_str());
    }
}

//! Get edge index for the vine (like 1, 2; 3)
//! @param e a descriptor for the edge.
//! @param tree a vine tree.
std::string StructureSelector::get_pc_index(
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

}
