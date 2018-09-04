
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

namespace vinecopulib {
    
// ------------------------- T-VINE STRUCTURE ---------------

TriangularArray<size_t> build_t_vine_array(RVineStructure struct0, size_t p)
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
        FamilySelector(data, selector.get_rvine_structure(), controls), 
        data_(data),
        d0_(vine_struct_.get_dim()),
        p_(d0_ / d_ - 1),
        lag_(0),
        struct0_(vine_struct_)
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

        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + d0_ * lag_;
            return index;
        };

        for (size_t t = 1; t < trees_.size(); t++) {
            auto old_tree = trees_[t];
            // add vertices for lagged variable
            for (auto v : boost::vertices(old_tree)) {
                auto v_new = boost::add_vertex(trees_[t]);
                
                // copy structure information
                trees_[t][v_new].conditioned = shift(trees_[t][v].conditioned);
                trees_[t][v_new].conditioning = shift(trees_[t][v].conditioning);
                trees_[t][v_new].all_indices = shift(trees_[t][v].all_indices);
                trees_[t][v_new].prev_edge_indices = shift(trees_[t][v].prev_edge_indices);
                
                // copy data and remove rows
                size_t n = trees_[t][v].hfunc1.rows() - 1;
                trees_[t][v_new].hfunc1 = trees_[t][v].hfunc1.bottomRows(n);
                trees_[t][v].hfunc1.conservativeResize(n);
                if (trees_[t][v].hfunc2.size() > 1) {
                    trees_[t][v_new].hfunc2 = trees_[t][v].hfunc2.bottomRows(n);
                    trees_[t][v].hfunc2.conservativeResize(n);
                }
            }
            
            // copy edges for lagged vine
            for (auto e : boost::edges(old_tree)) {
                size_t v1 = boost::source(e, old_tree);
                size_t v2 = boost::target(e, old_tree);
                auto e_new = boost::add_edge(v1 + d0_ * lag_, v2 + d0_ * lag_, trees_[t]).first;
                trees_[t][e_new].pair_copula = trees_[t][e].pair_copula;
                trees_[t][e_new].fit_id = trees_[t][e].fit_id;
            }
        }
        trees_opt_ = trees_;
        trees_ = std::vector<VineTree>(1);
        vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_),
                                      build_t_vine_array(struct0_, lag_));
        
        size_t n = data_.rows() - 1;
        Eigen::MatrixXd newdata(n, d_);
        newdata << data_.topRows(n), data_.rightCols(d0_).bottomRows(n);
        data_ = newdata;
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

    TVine(RVineStructure struct0, 
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
        vine_struct_ = RVineStructure(
            tools_stl::seq_int(1, d_),
            build_t_vine_array(struct0, p)
        );
        pair_copulas_ = make_pair_copula_store(d_);
    }
    
    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        if (vine_struct_.get_trunc_lvl() > 0) {
            auto newdata = data;
            auto revorder = vine_struct_.get_order();
            revorder.resize(d0_);
            tools_stl::reverse(revorder);
            for (size_t j = 0; j < d0_; ++j)
                newdata.col(j) = data.col(revorder[j] - 1);
            tools_select::TVineSelector selector(newdata, struct0_, controls);
            
            selector.select_all_trees(newdata);
            for (size_t lag = 1; lag <= p_; lag++) {
                selector.add_lag();
                selector.select_all_trees(selector.data());
            }
            
            vine_struct_ = selector.get_rvine_structure();
            threshold_ = selector.get_threshold();
            loglik_ = selector.get_loglik();
            nobs_ = data.rows();
            pair_copulas_ = selector.get_pair_copulas();
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
        
        vine_struct_ = tv_selector.get_rvine_structure();
        threshold_ = tv_selector.get_threshold();
        loglik_ = tv_selector.get_loglik();
        nobs_ = data.rows();
        pair_copulas_ = tv_selector.get_pair_copulas();
    }
    
private:
    void check_data_dim(const Eigen::MatrixXd &data) 
    { 
        if (d0_ != static_cast<size_t>(data.cols()))
            throw std::runtime_error("wrong number of columns");
    }
    
    size_t d0_;
    size_t p_;
    size_t in_;
    size_t out_;
    RVineStructure struct0_;
};



}
