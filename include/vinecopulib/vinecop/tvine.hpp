
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
    
// TODO: 
// - order selection
// - structure selection with fixed in or out node
// - loglikelihood calculation (current numbers are incorrect)
// - margins
// - truncated models
// - avoid unneccessary reordering
// 
// TO DISCUSS:
// - fit routine
// - reorder function
// ------------------------- T-VINE STRUCTURE ---------------

TriangularArray<size_t> get_actual_struct_array(const RVineStructure& structure)
{
    // structure array is in natural order, must convert to actual order
    auto no_array = structure.get_struct_array();
    auto order = structure.get_order();
    for (size_t i = 0; i < structure.get_dim() - 1; i++) {
        for (auto &s : no_array[i]) {
            s = order[s - 1];
        }
    }
    
    return no_array;
}

class TVineStructure : public RVineStructure {
public:    
    TVineStructure() {}
    
    TVineStructure(size_t cs_dim, size_t p = 0) : 
        TVineStructure(RVineStructure(tools_stl::seq_int(1, cs_dim)), p) 
    {}

    TVineStructure(const RVineStructure &cs_struct, 
                   size_t p = 0, 
                   size_t in_vertex = 1, 
                   size_t out_vertex = 1) :
        p_(p), 
        in_vertex_(in_vertex), 
        out_vertex_(out_vertex),
        cs_struct_(reorder_structure(cs_struct, in_vertex))
    {
        order_ = expand_order(cs_struct_.get_order(), p);
        struct_array_ = build_t_vine_array(cs_struct_, p, in_vertex, out_vertex);
        RVineStructure new_struct(order_, struct_array_);
        d_             = new_struct.get_dim();
        trunc_lvl_     = new_struct.get_trunc_lvl();
        struct_array_  = new_struct.get_struct_array();
        max_array_     = new_struct.get_max_array();
        needed_hfunc1_ = new_struct.get_needed_hfunc1();
        needed_hfunc2_ = new_struct.get_needed_hfunc2();    
    }
    
    size_t get_p() const
    {
        return p_;
    }
    
    size_t get_cs_dim() const
    {
        return cs_struct_.get_dim();
    }
        
    size_t get_in_vertex() const
    {
        return in_vertex_;
    }
    
    size_t get_out_vertex() const
    {
        return out_vertex_;
    }
    
    RVineStructure get_cs_structure() const 
    {
        return cs_struct_;
    }
    
    
private:
    std::vector<size_t> expand_order(const std::vector<size_t> &order, size_t p) const
    {
        size_t cs_dim = order.size();
        size_t d = cs_dim * (p + 1);
        std::vector<size_t> new_order(d);
        for (size_t i = 0; i < d; i++) {
            new_order[i] = order[i % cs_dim] + (i / cs_dim) * cs_dim;
        }
        
        return new_order;
    }
    
    std::vector<size_t> pivot_diag(const std::vector<size_t>& diag, 
                                   const TriangularArray<size_t>& struct_array, 
                                   size_t index)  const
    {
        size_t d = diag.size();
        std::vector<size_t> x(d);
        
        size_t i = 0;
        while (diag[i] != index)  {
            x[i] = diag[i];
            i++;
        }
        
        size_t pivot_col = i;
        while (i < d - 1) {
            x[i] = struct_array(d - 2 - i, pivot_col);
            i++;
        }
        x[d - 1] = index;
        
        return x;
    }
            
    RVineStructure reorder_structure(const RVineStructure& structure, size_t in_vertex) const
    {
        using namespace tools_stl;
        size_t d = structure.get_dim();
 
        // structure array is in natural order, must convert to actual order
        auto old_struct = get_actual_struct_array(structure);
        auto new_struct = old_struct;

        // prepare objects
        // working w/ rev_order is easier algorithmically; will be reversed again
        // when final object ist created.
        auto old_order = structure.get_rev_order();
        auto new_order = pivot_diag(old_order, old_struct, in_vertex);

        // loop through all columns
        for (size_t i = 0; i < d - 1; i++) {
            // extract elements of pivotal order that are required for the default
            // off-diagonal elements
            auto new_column = pivot_diag(old_order, old_struct, new_order[i]);
            new_column.erase(new_column.end() - 1);
            new_column.erase(new_column.begin(), new_column.begin() + i);
            tools_stl::reverse(new_column);
            
            auto diag_until = span(new_order, 0, i);
            auto diag_after = span(new_order, i, d - i);
            for (size_t t = 0; t < new_column.size(); t++) {
                // TODO: start from top
                // Check whether an element in this column is already contained in
                // the diagonal to the left. If so, we need to find another node 
                // that is connected to the diagonal entry of column i. We search 
                // for such an edge in the old structure, but only to the right of 
                // the in_vertex and to the the left of the current column (included).
                if (is_member(new_column[t], diag_until)) {
                    bool found_node = false;
                    // TODO: check <=
                    for (size_t j = find_position(in_vertex, old_order); j <= i; j++) {
                        if (new_order[i] == old_struct(t, j)) {
                            if (is_member(old_order[j], diag_after)) {
                                new_column[t] = old_order[j];
                                found_node = true;
                            }
                        } else if (new_order[i] == old_order[j]) {
                            if (is_member(old_struct(t, j), diag_after)) {
                                new_column[t] = old_struct(t, j);
                                found_node = true;
                            }
                        }
                        if (found_node) {
                            // The new entry may already be contained in this
                            // column. We need to check the next rows for that 
                            // as well.
                            diag_until = cat(new_column[t], diag_until);
                            break;
                        }
                    }
                }
            }
            new_struct[i] = new_column;
        }
        
        // this must always hold beacuse in_vertex comes last on the diagonal:
        new_struct(0, d - 2) = in_vertex;
        
        reverse(new_order);
        return RVineStructure(new_order, new_struct);
    }
   
    inline TriangularArray<size_t> build_t_vine_array(
        const RVineStructure& cs_struct, 
        size_t p,
        size_t in_vertex = 1,
        size_t out_vertex = 1) 
        const
    {
        size_t cs_dim = cs_struct.get_dim();
        size_t d = cs_dim * (p + 1);
        auto diag = cs_struct.get_rev_order();
        
        RVineStructure new_struct = cs_struct;
        if (diag[cs_dim - 1] != in_vertex) {
            new_struct = reorder_structure(new_struct, in_vertex);
        }
        
        auto struct_array = get_actual_struct_array(cs_struct);
 
        TriangularArray<size_t> strct(d);
        // copy cross-sectional structure
        for (size_t i = 0; i < cs_dim - 1; i++) {
            for (size_t j = 0; j < cs_dim - 1 - i; j++) {
                strct(i, j) = struct_array(i, j) + cs_dim * p;
            }
        }
        
        // fill parallelograms
        std::vector<size_t> par = pivot_diag(diag, struct_array, out_vertex);
        tools_stl::reverse(par);
        for (size_t lag = 1; lag <= p; lag++) {
            for (size_t i = 0; i < cs_dim; i++) {
                for (size_t j = 0; j < cs_dim; j++) {
                    strct(i + cs_dim * lag - j - 1, j) = par[i] + cs_dim * (p - lag);
                }
            }
        }

        // copy to other lags
        for (size_t lag = 1; lag <= p; lag++) {
            for (size_t j = 0; j < cs_dim; j++) {
                for (size_t i = 0; i < d - 1 - (j + cs_dim * lag); i++) {
                    strct(i, j + cs_dim * lag) = strct(i, j) - cs_dim * lag;
                }
            }
        }
        
        return strct;
    }
    
private:
    size_t p_;
    size_t in_vertex_;
    size_t out_vertex_;
    RVineStructure cs_struct_;
};


// ------------------------- SELECTOR ------------------------

inline Eigen::MatrixXd spread_lag(const Eigen::MatrixXd& data, size_t cs_dim)
{
    if (data.rows() < 2) {
        throw std::runtime_error("insufficient number of observations");
    }
    if (data.cols() % cs_dim != 0) {
        throw std::runtime_error("number of columns is not a multiple of cs_dim");
    }
    size_t n = data.rows() - 1;
    Eigen::MatrixXd newdata(n, data.cols() + cs_dim);
    newdata << data.topRows(n), data.rightCols(cs_dim).bottomRows(n);
    return newdata;
}

namespace tools_select {
    
class TVineSelector {
public:
    TVineSelector(const Eigen::MatrixXd &data, 
                  size_t in_vertex = 1, 
                  size_t out_vertex = 1) : 
        cs_dim_(data.cols()),
        in_vertex_(in_vertex),
        out_vertex_(out_vertex),
        data_(data)        
    {}

    size_t get_in_vertex() const
    {
        return in_vertex_;
    }
    
    size_t get_out_vertex() const
    {
        return out_vertex_;
    }
    
protected:
    size_t cs_dim_;
    size_t in_vertex_;
    size_t out_vertex_;
    Eigen::MatrixXd data_;
};
    
    
class TVineStructureSelector : public StructureSelector, public TVineSelector {
public:
    TVineStructureSelector(const Eigen::MatrixXd &data,
                           const FitControlsVinecop &controls) : 
        StructureSelector(data, controls), 
        TVineSelector(data)
    {}
        
    ~TVineStructureSelector() {}
        
    void select_tree(size_t t) override
    {
        auto new_tree = edges_as_vertices(trees_[t]);
        remove_edge_data(trees_[t]);
        add_allowed_edges(new_tree);
        if (boost::num_vertices(new_tree) > 2) {
            // has no effect in FamilySelector
            min_spanning_tree(new_tree);
        }
        if (boost::num_vertices(new_tree) > 0) {
                add_edge_info(new_tree);     // for pc estimation and next tree
                if (controls_.get_selection_criterion() == "mbicv") {
                    // adjust prior probability to tree level
                controls_.set_psi0(std::pow(psi0_, t + 1));
            }
            if (trees_opt_.size() > t + 1) {
                select_pair_copulas(new_tree, trees_opt_[t + 1]);
            } else {
                select_pair_copulas(new_tree);
            }
        }
        // make sure there is space for new tree
        trees_.resize(t + 2);
        trees_[t + 1] = new_tree;
    }
    
    void select_connecting_vertices()
    {
        double crit = 0.0, new_crit = 0.0;
        size_t d = data_.cols();
        size_t n = data_.rows() - 1;
        Eigen::MatrixXd pair_data(n, 2);
        for (size_t i = 0; i < d; i++) {
            for (size_t j = i; j < d; j++) {
                new_crit = wdm::wdm(data_.col(i).tail(n), 
                                    data_.col(j).head(n), 
                                    "hoeffd");
                if (std::abs(new_crit) > crit) {
                    crit = std::abs(new_crit);
                    in_vertex_ = i + 1;
                    out_vertex_ = j + 1;
                }
            }
        }
    }
    
protected:
    double compute_fit_id(const EdgeProperties& e) override
    {
        size_t min_c = std::min(e.conditioned[0], e.conditioned[1]);
        size_t max_c = std::max(e.conditioned[0], e.conditioned[1]);
        return (min_c % cs_dim_) * cs_dim_ * 10 + (max_c % cs_dim_);
    }
    
    void finalize(size_t trunc_lvl)
    {
        trees_opt_ = trees_;
        StructureSelector::finalize(trunc_lvl);
        trees_ = trees_opt_;
    }
};

class TVineFamilySelector : public FamilySelector, public TVineSelector {
public:
    TVineFamilySelector(const Eigen::MatrixXd &data,
                        const RVineStructure &cs_struct,
                        const FitControlsVinecop &controls,
                        size_t in_vertex = 1,
                        size_t out_vertex = 1)  : 
        FamilySelector(data, cs_struct, controls), 
        TVineSelector(data, in_vertex, out_vertex),
        lag_(0)
    {
        cs_struct_ = TVineStructure(cs_struct, 0, in_vertex, out_vertex);
    }
    
    TVineFamilySelector(const Eigen::MatrixXd &data,
                        const TVineStructureSelector &selector,
                        const FitControlsVinecop &controls) : 
        TVineFamilySelector(data, selector.get_rvine_structure(), controls,
                            selector.get_in_vertex(), selector.get_out_vertex())
    {
        trees_ = selector.get_trees();
        flip_pcs(trees_);
    }
    
    ~TVineFamilySelector() {}

    RVineStructure get_cs_structure() const
    {
        return cs_struct_;
    }
    
    void select_tree(size_t t) override
    {
        auto new_tree = edges_as_vertices(trees_[t]);
        remove_edge_data(trees_[t]);
        add_allowed_edges(new_tree);
        if (boost::num_vertices(new_tree) > 0) {
            add_edge_info(new_tree);       // for pc estimation and next tree
            if (controls_.get_selection_criterion() == "mbicv") {
                // adjust prior probability to tree level
                controls_.set_psi0(std::pow(psi0_, t + 1));
            }
            if (trees_opt_.size() > t + 1) {
                select_pair_copulas(new_tree, trees_opt_[t + 1]);
            } else {
                select_pair_copulas(new_tree);
            }
        }
        // make sure there is space for new tree
        trees_.resize(t + 2);
        trees_[t + 1] = new_tree;
    }
    
    void add_lag()
    {
        lag_++;
        d_ += cs_dim_;

        // add vertices and edges for lagged variable
        for (size_t t = 1; t < trees_.size(); t++) {
            auto old_tree = trees_[t];
            for (auto v : boost::vertices(old_tree)) 
                duplicate_vertex(v, trees_[t]);
            for (auto e : boost::edges(old_tree)) 
                duplicate_edge(e, trees_[t]);
        }

        // update trees and structure
        trees_opt_ = trees_;
        trees_ = std::vector<VineTree>(1);
        vine_struct_ = TVineStructure(cs_struct_, lag_, in_vertex_, out_vertex_);
        data_ = spread_lag(data_, cs_dim_);
    }
    
    Eigen::MatrixXd data()
    {
        return data_;
    }
    
protected:
    double compute_fit_id(const EdgeProperties& e) override
    {
        size_t min_c = std::min(e.conditioned[0], e.conditioned[1]);
        size_t max_c = std::max(e.conditioned[0], e.conditioned[1]);
        return (min_c % cs_dim_) * cs_dim_ * 10 + (max_c % cs_dim_);
    }
    
    void flip_pcs(std::vector<VineTree>& trees)
    {
        std::vector<size_t> rev_order = cs_struct_.get_rev_order();
        size_t d = rev_order.size();
        auto struct_array = get_actual_struct_array(cs_struct_);
        for (size_t t = 1; t < trees.size(); t++) {
            for (auto e : boost::edges(trees[t])) {
                std::vector<size_t> check_set(2);
                for (size_t k = 0; k < d; k++) {
                    // conditioned starts at 0 -> substract -1
                    check_set = {rev_order[k] - 1, struct_array(t - 1, k) - 1};
                    if (tools_stl::is_same_set(trees_[t][e].conditioned, check_set)) {
                        break;
                    }
                }
                if (trees[t][e].conditioned[0] != check_set[0]) {
                    trees[t][e].pair_copula.flip();
                }
            }
        }
    }
    
    void duplicate_vertex(size_t v, VineTree& tree)
    {
        auto v_new = boost::add_vertex(tree);
        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + cs_dim_ * lag_;
            return index;
        };

        // copy structure information
        tree[v_new].conditioned = shift(tree[v].conditioned);
        tree[v_new].conditioning = shift(tree[v].conditioning);
        tree[v_new].all_indices = shift(tree[v].all_indices);
        tree[v_new].prev_edge_indices = shift(tree[v].prev_edge_indices);
        
        // copy data and remove rows
        size_t n = tree[v].hfunc1.rows() - 1;
        tree[v_new].hfunc1 = tree[v].hfunc1.bottomRows(n);
        tree[v].hfunc1.conservativeResize(n);
        if (tree[v].hfunc2.size() > 1) {
            tree[v_new].hfunc2 = tree[v].hfunc2.bottomRows(n);
            tree[v].hfunc2.conservativeResize(n);
        }
    }
    
    void duplicate_edge(EdgeIterator e, VineTree& tree)
    {
        size_t v1 = boost::source(e, tree);
        size_t v2 = boost::target(e, tree);
        auto e_new = boost::add_edge(v1 + cs_dim_ * lag_, v2 + cs_dim_ * lag_, tree);
        tree[e_new.first].pair_copula = tree[e].pair_copula;
        tree[e_new.first].fit_id = tree[e].fit_id;
    }
    
    RVineStructure cs_struct_;
    size_t lag_;
};

} // end tools_select


// --------------------- T-VINE ----------------------------------
    
class TVine : public Vinecop {
public:
    TVine(size_t cs_dim, size_t p) : 
        TVine(RVineStructure(tools_stl::seq_int(1, cs_dim)), p) 
    {}

    TVine(const RVineStructure &cs_struct, 
          size_t p, 
          size_t in_vertex = 1, 
          size_t out_vertex = 1) : 
        cs_dim_(cs_struct.get_dim()),
        p_(p), 
        in_vertex_(in_vertex), 
        out_vertex_(out_vertex),
        tvine_struct_(TVineStructure(cs_struct, p, in_vertex, out_vertex))
    {
        threshold_ = 0.0;
        loglik_ = NAN;
	    d_ = tvine_struct_.get_dim();
        vine_struct_ = tvine_struct_;
        pair_copulas_ = make_pair_copula_store(d_);
    }
    
    TVine(const std::vector<std::vector<Bicop>> &pair_copulas,
          const TVineStructure &cs_struct,
          size_t p = 0,
          size_t in_vertex = 1,
          size_t out_vertex = 1) :
         TVine(cs_struct, p, in_vertex, out_vertex)
    {
        pair_copulas_ = pair_copulas;
    }
    
    size_t get_p() const
    {
        return p_;
    }
    
    size_t get_cs_dim() const
    {
        return cs_dim_;
    }
        
    size_t get_in_vertex() const
    {
        return in_vertex_;
    }
    
    size_t get_out_vertex() const
    {
        return out_vertex_;
    }
    
    RVineStructure get_cs_structure() const 
    {
        return tvine_struct_.get_cs_structure();
    }
    
    TVineStructure get_tvine_structure() const
    {
        return tvine_struct_;
    }
    
    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        if (vine_struct_.get_trunc_lvl() > 0) {
            tools_select::TVineFamilySelector selector(
                data, 
                tvine_struct_.get_cs_structure(), 
                controls,
                in_vertex_,
                out_vertex_);
            
            selector.select_all_trees(data);
            for (size_t lag = 1; lag <= p_; lag++) {
                selector.add_lag();
                selector.select_all_trees(selector.data());
            }
            
            finalize_fit(selector);
        }
    }
    
    void select_all(const Eigen::MatrixXd &data,
                    const FitControlsVinecop &controls)
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);
        
        tools_select::TVineStructureSelector selector(data, controls);
        selector.select_all_trees(data);
        selector.select_connecting_vertices();

        tools_select::TVineFamilySelector tv_selector(data, selector, controls);
        
        for (size_t lag = 1; lag <= p_; lag++) {
            tv_selector.add_lag();
            tv_selector.select_all_trees(tv_selector.data());
        }
        
        finalize_fit(tv_selector);
    }
    
    // Eigen::MatrixXd cond_simulate(size_t n, const Eigen::MatrixXd& data)
    // {
    //     check_data_dim(data);
    // 
    //     Eigen::MatrixXd cond_pit;
    //     Eigen::MatrixXd combined_vals(n, d_);
    //     if (p_ > 0) {
    //         // only most recent observations are used
    //         Eigen::MatrixXd cond_vals = data.bottomRows(p_);
    // 
    //         // spread data into columns
    //         for (size_t lag = 1; lag < p_; lag++) {
    //             cond_vals = spread_lag(cond_vals, cs_dim_);
    //         }
    // 
    //         d_ -= cs_dim_;
    //         vine_struct_.reduce(cond_vals.cols());
    //         combined_vals << 
    //             rosenblatt(cond_vals).replicate(n, 1), 
    //             tools_stats::simulate_uniform(n, cs_dim_);
    // 
    //         d_ += cs_dim_;
    //         vine_struct_ = TVine(cs_struct_, p_, in_, out_).get_rvine_structure();
    //     } else {
    //         combined_vals << tools_stats::simulate_uniform(n, cs_dim_);    
    //     }
    // 
    //     return inverse_rosenblatt(combined_vals).rightCols(cs_dim_);
    // }
    
    Eigen::MatrixXd simulate(const size_t n, 
                             const bool qrng = false, 
                             const size_t num_threads = 1,
                             const std::vector<int>& seeds = std::vector<int>()) const
    {
        Eigen::MatrixXd U;
        if (qrng) {
            if (cs_dim_ > 300) {
                U = tools_stats::sobol(n, cs_dim_, seeds);
            } else {
                U = tools_stats::ghalton(n, cs_dim_, seeds);
            }
        } else {
            U = tools_stats::simulate_uniform(n, cs_dim_, seeds);
        }
                
        // initialize first p + 1 lags
        Eigen::MatrixXd sim(n, cs_dim_);
        Eigen::MatrixXd Ui(1, d_);
        for (size_t i = 0; i <= p_; i++) {
            Ui.row(0).segment(i * cs_dim_, cs_dim_) = U.row(i);       
        }
        Eigen::MatrixXd V = inverse_rosenblatt(Ui, num_threads);
        for (size_t i = 0; i <= p_; i++) {
            sim.row(i) = V.block(0, i * cs_dim_, 1, cs_dim_);
        }
        
        // simulate conditional on previous observations
        for (size_t i = p_ + 1; i < n; i++) {
            Ui.leftCols(d_ - cs_dim_).swap(Ui.rightCols(d_ - cs_dim_));
            Ui.rightCols(cs_dim_) = U.row(i);
            sim.row(i) = inverse_rosenblatt(Ui, num_threads).rightCols(cs_dim_);
        }
        
        return sim;    
    }
    
    
protected:
    
    void finalize_tvine(const tools_select::TVineFamilySelector &selector)
    {
        in_vertex_ = selector.get_in_vertex();
        out_vertex_ = selector.get_out_vertex();
        tvine_struct_ = TVineStructure(selector.get_cs_structure(), 
                                       p_, in_vertex_, out_vertex_);
        Vinecop::finalize_fit(selector);
    }
    
    void check_data_dim(const Eigen::MatrixXd &data) const
    { 
        if (cs_dim_ != static_cast<size_t>(data.cols())) {
            std::stringstream msg;
            msg << "wrong number of columns." << std::endl <<
                "expected: " << cs_dim_ << std::endl <<
                "provided: " << data.cols() << std::endl;
            throw std::runtime_error(msg.str());
        }
    }
    
    size_t cs_dim_;
    size_t p_;
    size_t in_vertex_;
    size_t out_vertex_;
    TVineStructure tvine_struct_;
};

}
