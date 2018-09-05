
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
    
// ------------------------- T-VINE STRUCTURE ---------------

inline TriangularArray<size_t> build_t_vine_array(RVineStructure struct0, size_t p)
{
    size_t d0 = struct0.get_dim();
    size_t d = d0 * (p + 1);
    TriangularArray<size_t> strct(d);

    // copy cross-sectional structure
    for (size_t i = 0; i < d0 - 1; i++) {
        for (size_t j = 0; j < d0 - 1 - i; j++) {
            strct(i, j) = struct0.struct_array(i, j) + d0 * p;
        }
    }

    // fill parallelograms
    for (size_t lag = 1; lag <= p; lag++) {
        for (size_t i = 0; i < d0; i++) {
            for (size_t j = 0; j < d0; j++) {
                strct(i + d0 * lag - j - 1, j) = i + 1 + d0 * (p - lag);
            }
        }
    }

    // copy to other lags
    for (size_t lag = 1; lag <= p; lag++) {
        for (size_t j = 0; j < d0; j++) {
            for (size_t i = 0; i < d - 1 - (j + d0 * lag); i++) {
                strct(i, j + d0 * lag) = strct(i, j) - d0 * lag;
            }
        }
    }

    return strct;
}

inline Eigen::MatrixXd spread_lag(const Eigen::MatrixXd& data, size_t d0)
{
    if (data.rows() < 2) {
        throw std::runtime_error("insufficient number of observations");
    }
    if (data.cols() % d0 != 0) {
        throw std::runtime_error("number of columns is not a multiple of d0");
    }
    size_t n = data.rows() - 1;
    Eigen::MatrixXd newdata(n, data.cols() + d0);
    newdata << data.topRows(n), data.rightCols(d0).bottomRows(n);
    return newdata;
}

// ------------------------- SELECTOR ------------------------
namespace tools_select {

class TVineSelector : public FamilySelector {
public:
    TVineSelector(const Eigen::MatrixXd &data,
                  const RVineStructure &vine_struct,
                  const FitControlsVinecop &controls)  : 
        FamilySelector(data, vine_struct, controls), 
        data_(data),
        d0_(vine_struct.get_dim()), 
        p_(d0_ / d_ - 1), 
        lag_(0), 
        struct0_(vine_struct)
    {}
    
    TVineSelector(const Eigen::MatrixXd &data,
                  const VinecopSelector &selector,
                  const FitControlsVinecop &controls) : 
        TVineSelector(data, selector.get_rvine_structure(), controls)
    {
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
        d_ += d0_;

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
        vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_),
                                      build_t_vine_array(struct0_, lag_));
        data_ = spread_lag(data_, d0_);
    }
    
    void duplicate_vertex(size_t v, VineTree& tree)
    {
        auto v_new = boost::add_vertex(tree);
        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + d0_ * lag_;
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
        auto e_new = boost::add_edge(v1 + d0_ * lag_, v2 + d0_ * lag_, tree);
        tree[e_new.first].pair_copula = tree[e].pair_copula;
        tree[e_new.first].fit_id = tree[e].fit_id;
    }
        
    Eigen::MatrixXd data()
    {
        return data_;
    }
    
protected:
    double compute_fit_id(const EdgeProperties& e) override
    {
        return (e.conditioned[0] % d0_) * d0_ * 10 + (e.conditioned[1] % d0_);
    }
    
    Eigen::MatrixXd data_;
    size_t d0_;
    size_t p_;
    size_t lag_;
    RVineStructure struct0_;
};

} // end tools_select


// --------------------- T-VINE ----------------------------------
    
class TVine : public Vinecop {
public:
    TVine(size_t d0, size_t p = 0) : 
        TVine(RVineStructure(tools_stl::seq_int(1, d0)), p) {}

    TVine(const RVineStructure &struct0, 
          size_t p = 0, 
          size_t in = 1, 
          size_t out = 1) : 
        Vinecop(struct0), 
        d0_(struct0.get_dim()),
        p_(p), 
        in_(in), 
        out_(out), 
        struct0_(struct0)
    {
		d_ = struct0.get_dim() * (p + 1);
        std::vector<size_t> order(d_);
        for (size_t i = 0; i < d_; i++) {
            order[i] = struct0_.get_order()[i % d0_] + (i / d0_) * d0_;
        }
        vine_struct_ = RVineStructure(order, build_t_vine_array(struct0, p));
        pair_copulas_ = make_pair_copula_store(d_);
    }
    
    TVine(const std::vector<std::vector<Bicop>> &pair_copulas,
          const RVineStructure &vine_struct, 
          size_t p) : 
         Vinecop(vine_struct), 
         d0_(vine_struct.get_dim() / (p + 1)),
         p_(p), 
         in_(1), 
         out_(1)
    {
          struct0_ = vine_struct;
          struct0_.reduce(d0_);
          check_pair_copulas_rvine_structure(pair_copulas);
          pair_copulas_ = pair_copulas;
    }
    
    RVineStructure get_struct0() const
    {
        return struct0_;
    }
    
    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        if (vine_struct_.get_trunc_lvl() > 0) {
            auto newdata = data;
            auto rev_order = vine_struct_.get_order();
            rev_order.resize(d0_);
            tools_stl::reverse(rev_order);
            for (size_t j = 0; j < d0_; ++j)
                newdata.col(j) = data.col(rev_order[j] - 1);
            tools_select::TVineSelector selector(newdata, struct0_, controls);
            
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
    
    Eigen::MatrixXd cond_simulate(size_t n, const Eigen::MatrixXd& data)
    {
        check_data_dim(data);
        
        Eigen::MatrixXd cond_pit;
        Eigen::MatrixXd combined_vals(n, d_);
        if (p_ > 0) {
            // only most recent observations are used
            Eigen::MatrixXd cond_vals = data.bottomRows(p_);
            
            // spread data into columns
            for (size_t lag = 1; lag < p_; lag++) {
                cond_vals = spread_lag(cond_vals, d0_);
            }
            
            d_ -= d0_;
            vine_struct_.reduce(cond_vals.cols());
            combined_vals << 
                rosenblatt(cond_vals).replicate(n, 1), 
                tools_stats::simulate_uniform(n, d0_);
            
            d_ += d0_;
            vine_struct_ = TVine(struct0_, p_, in_, out_).get_rvine_structure();
        } else {
            combined_vals << tools_stats::simulate_uniform(n, d0_);    
        }
        
        return inverse_rosenblatt(combined_vals).rightCols(d0_);
    }
    
    Eigen::MatrixXd simulate(const size_t n, 
                             const bool qrng = false, 
                             const size_t num_threads = 1,
                             const std::vector<int>& seeds = std::vector<int>()) const
    {
        auto uniform_rng = [seeds, qrng] (size_t n, size_t d) {
            if (qrng) {
                if (d > 300) {
                    return tools_stats::sobol(n, d, seeds);
                } else {
                    return tools_stats::ghalton(n, d, seeds);
                }
            } else {
                return tools_stats::simulate_uniform(n, d, seeds);
            }
        };
        
        // initialize
        Eigen::MatrixXd sim(n, d0_);
        Eigen::MatrixXd Ui = uniform_rng(1, d_);
        Eigen::MatrixXd V = inverse_rosenblatt(Ui, num_threads);
        for (size_t i = 0; i <= p_; i++) {
            sim.row(i) = V.block(0, i * d0_, 1, d0_);
        }
        
        // simulate conditional on previous observations
        Eigen::MatrixXd U = uniform_rng(n, d0_);
        for (size_t i = p_ + 1; i < n; i++) {
            Ui.leftCols(d_ - d0_).swap(Ui.rightCols(d_ - d0_));
            Ui.rightCols(d0_) = U.row(i);
            sim.row(i) = inverse_rosenblatt(Ui, num_threads).rightCols(d0_);
        }
        
        return sim;    
    }
    
    
private:
    void check_data_dim(const Eigen::MatrixXd &data) 
    { 
        if (d0_ != static_cast<size_t>(data.cols())) {
            std::stringstream msg;
            msg << "cond_vals has wrong number of columns." << std::endl <<
                "expected: " << d0_ << std::endl <<
                "provided: " << data.cols() << std::endl;
            throw std::runtime_error(msg.str());
        }
    }
    
    size_t d0_;
    size_t p_;
    size_t in_;
    size_t out_;
    RVineStructure struct0_;
};

}
