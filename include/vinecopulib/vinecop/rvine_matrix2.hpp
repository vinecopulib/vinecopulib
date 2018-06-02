#pragma once

#include <Eigen/Dense>
#include <vector>
#include <chrono>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

template<class T>
class RVineMatrix2 {
public:
    RVineMatrix2() {}
    RVineMatrix2(size_t d) : d_(d)
    {
        mat_ = std::vector<std::vector<T>>(d - 1);
        for (size_t i = 0; i < mat_.size(); i++)
            mat_[i] = std::vector<T>(mat_.size() - i, 0);
    }
    
    T& operator()(size_t tree, size_t edge) {return mat_[edge][tree];}
    T operator()(size_t tree, size_t edge) const {return mat_[edge][tree];}

    std::string str()
    {
        std::stringstream str;
        for (size_t i = 0; i < mat_.size(); i++) {
            for (size_t j = 0; j < mat_.size() - i; j++) {
                str << (*this)(i, j) << " ";
            }
            str << std::endl;
        }
        return str.str();
    }
    
    size_t dim() const {return d_;}

private:
    size_t d_;
    std::vector< std::vector<T> > mat_;
};


class RVineStructure {
public:
    RVineStructure() {}
    
    RVineStructure(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat);

    std::vector<size_t> get_order() const {return order_;}
    RVineMatrix2<size_t> get_struct_matrix() const {return mat_;}
    RVineMatrix2<size_t> get_max_matrix() const {return max_mat_;}
    RVineMatrix2<size_t> get_needed_hfunc1() const {return needed_hfunc1_;}
    RVineMatrix2<size_t> get_needed_hfunc2() const {return needed_hfunc2_;}

    size_t struct_matrix(size_t tree, size_t edge) const {
        return mat_(tree, edge);
    }

    size_t max_matrix(size_t tree, size_t edge) const {
        return max_mat_(tree, edge);
    }

    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const {
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix(d_, d_);
        matrix.fill(0);
        for (size_t i = 0; i < d_; ++i) {
            for (size_t j = 0; j < d_ - i - 1; ++j) {
                matrix(i, j) = order_[mat_(i, j) - 1];
            }
            matrix(d_ - i - 1, i) = order_[d_ - i - 1];
        }
        return matrix;
    }

protected:

    size_t find_trunc_lvl(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat)
    {
        size_t trunc_lvl; ;
        for (trunc_lvl = mat.cols() - 1; trunc_lvl > 0; trunc_lvl--) {
            if (mat(trunc_lvl - 1, 0) != 0)
                break;
        }

        return trunc_lvl;
    }

    std::vector<size_t> compute_order(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
    {
        std::vector<size_t> order(d_);
        for (size_t i = 0; i < d_; i++)
            order[i] = mat(i, d_ - i - 1);
        
        return order;
    }
    
    RVineMatrix2<size_t> compute_struct_mat(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
    {

        // create vector of new variable labels
        auto order = tools_stl::get_order(get_order());

        // copy upper triangle and relabel to natural order
        RVineMatrix2<size_t> struct_mat(d_);
        for (size_t i = 0; i < d_ - 1; i++) {
            for (size_t j = 0; j < d_ - 1 - i; j++) {
                struct_mat(i, j) = order[mat(i, j) - 1] + 1;
            }
        }

        return struct_mat;
    }
    
    RVineMatrix2<size_t> compute_max_mat() const
    {
        RVineMatrix2<size_t> max_mat = mat_;
        for (size_t j = 0; j < d_ - 1; j++) {
            for (size_t i = 1; i < d_ - 1 - j; i++) {
                    max_mat(i, j) = std::max(mat_(i, j), max_mat(i - 1, j));
            }
        }
    
        return max_mat;
    }
    
    RVineMatrix2<size_t> compute_needed_hfunc1() const
    {
        RVineMatrix2<size_t> needed_hfunc1(d_);
 
        for (size_t i = 0; i < d_ - 2; i++) {
            for (size_t j = 0; j < d_ - 2 - i; j++) {
                if (mat_(i + 1, j) != max_mat_(i + 1, j))
                    needed_hfunc1(i, d_ - max_mat_(i + 1, j)) = 1;
            }
        }
        
        return needed_hfunc1;
    }
    
    RVineMatrix2<size_t> compute_needed_hfunc2() const
    {
        RVineMatrix2<size_t> needed_hfunc2(d_);
 
        for (size_t i = 0; i < d_ - 2; i++) {
            for (size_t j = 0; j < d_ - 2 - i; j++) {
                needed_hfunc2(i, j) = 1;
                if (mat_(i + 1, j) == max_mat_(i + 1, j))
                    needed_hfunc2(i, d_ - max_mat_(i + 1, j)) = 1;
            }
        }
                
        return needed_hfunc2;
    }


private:
    std::vector<size_t> order_;
    size_t d_;
    size_t trunc_lvl_;
    RVineMatrix2<size_t> mat_;
    RVineMatrix2<size_t> max_mat_;
    RVineMatrix2<size_t> needed_hfunc1_;
    RVineMatrix2<size_t> needed_hfunc2_;
};

inline RVineStructure::RVineStructure(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
    d_ = mat.cols();
    trunc_lvl_ = find_trunc_lvl(mat);
    order_ = compute_order(mat);
    mat_ = compute_struct_mat(mat);
    max_mat_ = compute_max_mat();
    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
}

}

class Timer {
public:
    void start()
    {
        start_ = std::chrono::steady_clock::now();
    }
    
    void end() 
    {
        end_ = std::chrono::steady_clock::now();
        auto diff = end_ - start_;
        std::cout << diff.count() << std::endl;
    }
    
    
private:
    std::chrono::time_point<std::chrono::steady_clock> start_;
    std::chrono::time_point<std::chrono::steady_clock> end_;
};
