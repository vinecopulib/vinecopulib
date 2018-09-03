
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

namespace vinecopulib {
    
// ------------------------- SHIFT OPERATOR ---------------------------
class LagShift {
public:
    LagShift(size_t d, size_t lag) : d(d), lag(lag) {}
    
    friend size_t operator+(const size_t& index, const LagShift& shift); 
    friend size_t operator-(const size_t& index, const LagShift& shift); 
    friend std::vector<size_t>  operator+(const std::vector<size_t>& index, 
                                          const LagShift& shift); 
    friend std::vector<size_t>  operator-(const std::vector<size_t>& index, 
                                          const LagShift& shift); 

    size_t d;
    size_t lag;
};

size_t operator+(const size_t& index, const LagShift& shift)
{
    return index + shift.d * shift.lag;
}
size_t operator-(const size_t& index, const LagShift& shift)
{
    assert(index >= shift.d * shift.lag);
    return index - shift.d * shift.lag;
}

std::vector<size_t> operator+(const std::vector<size_t>& index, const LagShift& shift)
{
    auto new_index = index;
    for (auto& i : new_index)
        i = i + shift;
    return new_index;
}

std::vector<size_t> operator-(const std::vector<size_t>& index, const LagShift& shift)
{
    auto new_index = index;
    for (auto& i : new_index)
        i = i - shift;
    return new_index;
}


// ------------------------- T-VINE STRUCTURE ---------------

TriangularArray<size_t> build_t_vine_array(RVineStructure cs_struct_, size_t order_)
{
    size_t cs_dim_ = cs_struct_.get_dim();
    size_t d_ = cs_dim_ * order_;
    TriangularArray<size_t> strct(d_);
    LagShift val_shift(cs_dim_, order_ - 1);

    // copy cross-sectional structure
    for (size_t i = 0; i < cs_dim_ - 1; i++) {
        for (size_t j = 0; j < cs_dim_ - 1 - i; j++) {
            strct(i, j) = cs_struct_.struct_array(i, j) + val_shift;
        }
    }

    // fill parallelograms
    for (size_t lag = 1; lag < order_; lag++) {
        val_shift = LagShift(cs_dim_, order_ - lag - 1);
        LagShift shift(cs_dim_, lag);
        for (size_t i = 0; i < cs_dim_; i++) {
            for (size_t j = 0; j < cs_dim_; j++) {
                strct((i + shift) - j - 1, j) = (i + 1) + val_shift;
            }
        }
    }

    // copy to other lags
    for (size_t lag = 1; lag < order_; lag++) {
        LagShift shift(cs_dim_, lag);
        for (size_t j = 0; j < cs_dim_; j++) {
            for (size_t i = 0; i < d_ - 1 - (j + shift); i++) {
                strct(i, j + shift) = strct(i, j) - shift;
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
                  const FitControlsVinecop &controls) 
                  : 
                  FamilySelector(data, vine_struct, controls),
                  order_(vine_struct.get_dim() / d_),
                  lag_(0),
                  cs_struct_(vine_struct),
                  cs_dim_(vine_struct.get_dim())
    {}
    
    TVineSelector(const Eigen::MatrixXd &data,
                  const VinecopSelector &selector,
                  const FitControlsVinecop &controls) 
                  : 
                  cs_struct_(selector.get_rvine_structure()),
                  FamilySelector(data, selector.get_rvine_structure(), controls),
                  order_(selector.get_rvine_structure().get_dim() / d_),
                  lag_(0),
                  cs_dim_(selector.get_rvine_structure().get_dim())
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
        LagShift lag(d_, ++lag_);
        for (size_t t = 1; t < trees_.size(); t++) {
            auto old_tree = trees_[t];
            // add vertices for lagged variable
            for (auto v : boost::vertices(old_tree)) {
                auto v_new = boost::add_vertex(trees_[t]);
                
                // copy structure information
                trees_[t][v_new].conditioned = trees_[t][v].conditioned + lag;
                trees_[t][v_new].conditioning = trees_[t][v].conditioning + lag;
                trees_[t][v_new].all_indices = trees_[t][v].all_indices + lag;
                trees_[t][v_new].prev_edge_indices = trees_[t][v].prev_edge_indices + lag;
                
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
                auto e_new = boost::add_edge(v1 + lag, v2 + lag, trees_[t]).first;
                trees_[t][e_new].pair_copula = trees_[t][e].pair_copula;
                trees_[t][e_new].fit_id = trees_[t][e].fit_id;
            }
        }
        trees_opt_ = trees_;
        trees_ = std::vector<VineTree>(1);
        d_ += cs_dim_; 
        vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_),
                                      build_t_vine_array(cs_struct_, (lag_ + 1)));
    }
    
    Eigen::MatrixXd add_lag(const Eigen::MatrixXd& data)
    {
        size_t n = data.rows() - 1;
        Eigen::MatrixXd newdata(n, d_);
        newdata << data.topRows(n), data.rightCols(cs_dim_).bottomRows(n);
        return newdata;
    }
    
protected:
    double compute_fit_id(const EdgeProperties& e) override
    {
        return (e.conditioned[0] % cs_dim_) * cs_dim_ * 10 + (e.conditioned[1] % cs_dim_);
    }
    
    size_t order_;
    size_t lag_;
    RVineStructure cs_struct_;
    size_t cs_dim_;
};

} // end tools_select
    
class TVine : public Vinecop {
public:
    TVine(RVineStructure cs_struct, 
          size_t order = 1, 
          size_t in = 1, 
          size_t out = 1) 
          : 
          cs_dim_(cs_struct.get_dim()),
          order_(order),
          in_(in),
          out_(out),
          cs_struct_(cs_struct),
          Vinecop(cs_struct)
    {
		d_ = cs_struct.get_dim() * order;
        vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_), 
                                      build_t_vine_array(cs_struct, order));
        pair_copulas_ = make_pair_copula_store(d_);
    }
    
    TVine(size_t d, size_t order = 1) : 
        TVine(RVineStructure(tools_stl::seq_int(1, d)), order) {}
    
    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        if (vine_struct_.get_trunc_lvl() > 0) {
            auto newdata = data;
            auto revorder = vine_struct_.get_order();
            revorder.resize(cs_dim_);
            tools_stl::reverse(revorder);
            for (size_t j = 0; j < cs_dim_; ++j)
                newdata.col(j) = data.col(revorder[j] - 1);
            tools_select::TVineSelector selector(newdata, cs_struct_, controls);
            
            selector.select_all_trees(newdata);
            for (size_t lag = 1; lag < order_; lag++) {
                selector.add_lag();
                newdata = selector.add_lag(newdata);
                selector.select_all_trees(newdata);
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
        for (size_t lag = 1; lag < order_; lag++) {
            tv_selector.add_lag();
            newdata = tv_selector.add_lag(newdata);
            tv_selector.select_all_trees(newdata);
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
        if (cs_dim_ != data.cols())
            throw std::runtime_error("wrong number of columns");
    }
        

    size_t cs_dim_;
    size_t order_;
    size_t in_;
    size_t out_;
    RVineStructure cs_struct_;
};



}
