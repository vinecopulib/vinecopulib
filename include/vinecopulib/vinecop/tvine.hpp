
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
// 
// TO DISCUSS:
// - fit routine
// - reorder function
// ------------------------- T-VINE STRUCTURE ---------------

class TVineStructure : public RVineStructure {
public:
    TVineStructure() {}
    
    TVineStructure(size_t d_cs, size_t p = 0) : 
        TVineStructure(RVineStructure(tools_stl::seq_int(1, d_cs)), p) 
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
        d_ = order_.size();
        struct_array_ = build_t_vine_array(cs_struct_, p, in_vertex, out_vertex);
        to_natural_order();
        trunc_lvl_ = struct_array_.get_trunc_lvl();
        max_array_     = compute_max_array();
        needed_hfunc1_ = compute_needed_hfunc1();
        needed_hfunc2_ = compute_needed_hfunc2();
    }
    
    size_t get_p() const
    {
        return p_;
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
        size_t d_cs = order.size();
        size_t d = d_cs * (p + 1);
        std::vector<size_t> new_order(d);
        for (size_t i = 0; i < d; i++) {
            new_order[i] = order[i % d_cs] + (i / d_cs) * d_cs;
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
    
    TriangularArray<size_t> get_actual_struct_array(const RVineStructure& structure) const
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
                // Check whether an element in this column is already contained in
                // the diagonal to the left. If so, we need to find another node 
                // that is connected to the diagonal entry of column i. We search 
                // for such an edge in the old structure, but only to the right of 
                // the in_vertex and to the the left of the current column (included).
                if (is_member(new_column[t], diag_until)) {
                    bool found_node = false;
                    for (size_t j = find_position(in_vertex, old_order); j <= i; j++) {
                        if (new_order[i] == old_struct(t, j)) {
                            if (is_member(old_order[j], diag_after)) {
                                new_column[t] = old_order[j];
                                found_node = true;
                            }
                        } else if (new_order[i] == old_order[j]) {
                            if (is_member(old_struct(t, j), diag_after)) {
                                new_column[t] = old_struct(t, j);
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
        size_t d_cs = cs_struct.get_dim();
        size_t d = d_cs * (p + 1);
        auto diag = cs_struct.get_rev_order();
        
        RVineStructure new_struct = cs_struct;
        if (diag[d_cs - 1] != in_vertex) {
            new_struct = reorder_structure(new_struct, in_vertex);
        }
        
        auto struct_array = get_actual_struct_array(cs_struct);
 
        TriangularArray<size_t> strct(d);
        // copy cross-sectional structure
        for (size_t i = 0; i < d_cs - 1; i++) {
            for (size_t j = 0; j < d_cs - 1 - i; j++) {
                strct(i, j) = struct_array(i, j) + d_cs * p;
            }
        }
        
        // fill parallelograms
        std::vector<size_t> par = pivot_diag(diag, struct_array, out_vertex);
        tools_stl::reverse(par);
        for (size_t lag = 1; lag <= p; lag++) {
            for (size_t i = 0; i < d_cs; i++) {
                for (size_t j = 0; j < d_cs; j++) {
                    strct(i + d_cs * lag - j - 1, j) = par[i] + d_cs * (p - lag);
                }
            }
        }

        // copy to other lags
        for (size_t lag = 1; lag <= p; lag++) {
            for (size_t j = 0; j < d_cs; j++) {
                for (size_t i = 0; i < d - 1 - (j + d_cs * lag); i++) {
                    strct(i, j + d_cs * lag) = strct(i, j) - d_cs * lag;
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

inline Eigen::MatrixXd spread_lag(const Eigen::MatrixXd& data, size_t d_cs)
{
    if (data.rows() < 2) {
        throw std::runtime_error("insufficient number of observations");
    }
    if (data.cols() % d_cs != 0) {
        throw std::runtime_error("number of columns is not a multiple of d_cs");
    }
    size_t n = data.rows() - 1;
    Eigen::MatrixXd newdata(n, data.cols() + d_cs);
    newdata << data.topRows(n), data.rightCols(d_cs).bottomRows(n);
    return newdata;
}

namespace tools_select {

class TVineSelector : public FamilySelector {
public:
    TVineSelector(const Eigen::MatrixXd &data,
                  const RVineStructure &cs_struct,
                  const FitControlsVinecop &controls,
                  size_t in_vertex = 1,
                  size_t out_vertex = 1)  : 
        FamilySelector(data, cs_struct, controls), 
        data_(data),
        d_cs_(cs_struct.get_dim()), 
        p_(d_cs_ / d_ - 1),
        in_vertex_(in_vertex),
        out_vertex_(out_vertex),
        lag_(0),
        cs_struct_(cs_struct)
    {}
    
    TVineSelector(const Eigen::MatrixXd &data,
                  const VinecopSelector &selector,
                  const FitControlsVinecop &controls) : 
        TVineSelector(data, selector.get_rvine_structure(), controls)
    {
        select_connecting_vertices();
        pair_copulas_ = selector.get_pair_copulas();
        trees_ = selector.get_trees_opt();      
    }
    
    ~TVineSelector() {}
    
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
        d_ += d_cs_;

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
        vine_struct_ = TVineStructure(cs_struct_, lag_);
        data_ = spread_lag(data_, d_cs_);
    }
    
    Eigen::MatrixXd data()
    {
        return data_;
    }
    
protected:
    double compute_fit_id(const EdgeProperties& e) override
    {
        return (e.conditioned[0] % d_cs_) * d_cs_ * 10 + (e.conditioned[1] % d_cs_);
    }
    
private:
    void select_connecting_vertices()
    {
        double crit = 0.0, new_crit = 0.0;
        size_t d = data_.cols();
        size_t n = data_.rows() - 1;
        Eigen::MatrixXd pair_data(n, 2);
        for (size_t i = 0; i < d; i++) {
            for (size_t j = i; j < d; j++) {
                new_crit = wdm::wdm(data_.col(i).tail(n), data_.col(j).head(n), "hoeffd");
                if (std::abs(new_crit) > crit) {
                    crit = std::abs(new_crit);
                    in_vertex_ = i;
                    out_vertex_ = j;
                }
            }
        }
    }
    
    void duplicate_vertex(size_t v, VineTree& tree)
    {
        auto v_new = boost::add_vertex(tree);
        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + d_cs_ * lag_;
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
        auto e_new = boost::add_edge(v1 + d_cs_ * lag_, v2 + d_cs_ * lag_, tree);
        tree[e_new.first].pair_copula = tree[e].pair_copula;
        tree[e_new.first].fit_id = tree[e].fit_id;
    }

    
    Eigen::MatrixXd data_;
    size_t d_cs_;
    size_t p_;
    size_t in_vertex_;
    size_t out_vertex_;
    size_t lag_;
    RVineStructure cs_struct_;
};

} // end tools_select


// --------------------- T-VINE ----------------------------------
    
class TVine : public Vinecop {
public:
    TVine(size_t d_cs, size_t p = 0) : 
        TVine(RVineStructure(tools_stl::seq_int(1, d_cs)), p) {}

    TVine(const RVineStructure &cs_struct, 
          size_t p = 0, 
          size_t in_vertex = 1, 
          size_t out_vertex = 1) : 
        d_cs_(cs_struct.get_dim()),
        p_(p), 
        in_vertex_(in_vertex), 
        out_vertex_(out_vertex)
    {
        threshold_ = 0.0;
        loglik_ = NAN;
	    d_ = d_cs_ * (p + 1);
        vine_struct_ = TVineStructure(cs_struct, p, in_vertex, out_vertex);
        pair_copulas_ = make_pair_copula_store(d_);
    }
    
    TVine(const std::vector<std::vector<Bicop>> &pair_copulas,
          const TVineStructure &tvine_struct, 
          size_t p) : 
         d_cs_(tvine_struct.get_dim() / (p + 1)),
         p_(p), 
         in_vertex_(1), 
         out_vertex_(1)
    {
        d_ = tvine_struct.get_dim();
        threshold_ = 0.0;
        loglik_ = NAN;
        check_pair_copulas_rvine_structure(pair_copulas);
        pair_copulas_ = pair_copulas;
    }
    
    size_t get_p() const
    {
        return p_;
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
        return vine_struct_.get_cs_structure();
    }
    
    TVineStructure get_rvine_structure() const
    {
        return vine_struct_;
    }
    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        if (vine_struct_.get_trunc_lvl() > 0) {
            auto newdata = data;
            auto rev_order = vine_struct_.get_order();
            rev_order.resize(d_cs_);
            tools_stl::reverse(rev_order);
            for (size_t j = 0; j < d_cs_; ++j)
                newdata.col(j) = data.col(rev_order[j] - 1);
            tools_select::TVineSelector selector(newdata, 
                                                 vine_struct_.get_cs_structure(), 
                                                 controls,
                                                 in_vertex_, 
                                                 out_vertex_);
            
            selector.select_all_trees(newdata);
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
        
        tools_select::StructureSelector selector(data, controls);
        selector.select_all_trees(data);

        auto newdata = data;
        tools_select::TVineSelector tv_selector(data, selector, controls);
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
    //             cond_vals = spread_lag(cond_vals, d_cs_);
    //         }
    // 
    //         d_ -= d_cs_;
    //         vine_struct_.reduce(cond_vals.cols());
    //         combined_vals << 
    //             rosenblatt(cond_vals).replicate(n, 1), 
    //             tools_stats::simulate_uniform(n, d_cs_);
    // 
    //         d_ += d_cs_;
    //         vine_struct_ = TVine(cs_struct_, p_, in_, out_).get_rvine_structure();
    //     } else {
    //         combined_vals << tools_stats::simulate_uniform(n, d_cs_);    
    //     }
    // 
    //     return inverse_rosenblatt(combined_vals).rightCols(d_cs_);
    // }
    
    Eigen::MatrixXd simulate(const size_t n, 
                             const bool qrng = false, 
                             const size_t num_threads = 1,
                             const std::vector<int>& seeds = std::vector<int>()) const
    {
        Eigen::MatrixXd U;
        if (qrng) {
            if (d_cs_ > 300) {
                U = tools_stats::sobol(n, d_cs_, seeds);
            } else {
                U = tools_stats::ghalton(n, d_cs_, seeds);
            }
        } else {
            U = tools_stats::simulate_uniform(n, d_cs_, seeds);
        }
        
        // initialize first p + 1 lags
        Eigen::MatrixXd sim(n, d_cs_);
        Eigen::MatrixXd Ui(1, d_);
        for (size_t i = 0; i <= p_; i++) {
            Ui << U.row(i);       
        }
        Eigen::MatrixXd V = inverse_rosenblatt(Ui, num_threads);
        for (size_t i = 0; i <= p_; i++) {
            sim.row(i) = V.block(0, i * d_cs_, 1, d_cs_);
        }
        
        // simulate conditional on previous observations
        for (size_t i = p_ + 1; i < n; i++) {
            Ui.leftCols(d_ - d_cs_).swap(Ui.rightCols(d_ - d_cs_));
            Ui.rightCols(d_cs_) = U.row(i);
            sim.row(i) = inverse_rosenblatt(Ui, num_threads).rightCols(d_cs_);
        }
        
        return sim;    
    }
    
    
private:

    void check_data_dim(const Eigen::MatrixXd &data) const
    { 
        if (d_cs_ != static_cast<size_t>(data.cols())) {
            std::stringstream msg;
            msg << "wrong number of columns." << std::endl <<
                "expected: " << d_cs_ << std::endl <<
                "provided: " << data.cols() << std::endl;
            throw std::runtime_error(msg.str());
        }
    }
    
    size_t d_cs_;
    size_t p_;
    size_t in_vertex_;
    size_t out_vertex_;
    TVineStructure vine_struct_;
};

}
