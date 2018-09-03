
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

namespace vinecopulib {
    
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
          out_(out)
    {
		d_ = cs_struct.get_dim() * order;
        TriangularArray<size_t> strct(d_);
        
        // copy cross-sectional structure
        for (size_t i = 0; i < cs_dim_ - 1; i++) {
            for (size_t j = 0; j < cs_dim_ - 1 - i; j++) {
             	strct(i, j) = shift(cs_struct.struct_array(i, j), order - 1);
            }
        }

        // fill parallelograms
        for (size_t lag = 1; lag < order; lag++) {
            for (size_t i = 0; i < cs_dim_; i++) {
                for (size_t j = 0; j < cs_dim_; j++) {
                    strct(shift(i, lag) - j - 1, j) = shift(i + 1, order - lag - 1);
                }
            }
        }

        // copy to other lags
        for (size_t lag = 1; lag < order; lag++) {
        	for (size_t j = 0; j < cs_dim_; j++) {
                  for (size_t i = 0; i < d_ - 1 - shift(j, lag); i++) {
                    strct(i, shift(j, lag)) =  shift_back(strct(i, j), lag);
                }
            }
        }

        vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_), strct);
        pair_copulas_ = make_pair_copula_store(d_);
    }
    
    
private:
    size_t shift(size_t index, size_t lag) {return index + cs_dim_ * lag;}
    size_t shift_back(size_t index, size_t lag) {return index - cs_dim_ * lag;}

    size_t cs_dim_;
    size_t order_;
    size_t in_;
    size_t out_;
    RVineStructure vine_struct_;
};

namespace tools_select {

class TVineSelector : public VinecopSelector {
public:
    TVineSelector(const Eigen::MatrixXd &data,
                  const RVineStructure &vine_struct,
                  const FitControlsVinecop &controls);
    ~TVineSelector()
    {
    }

private:
    void add_allowed_edges(VineTree &tree);
    void finalize(size_t trunc_lvl);
};

}

}
