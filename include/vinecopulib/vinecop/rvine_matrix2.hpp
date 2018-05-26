#pragma once

#include <Eigen/Dense>
#include <vector>
#include <chrono>

template<class T>
class RVineMatrix2 {
public:
    RVineMatrix2() {}
    RVineMatrix2(size_t d) : d_(d)
    {
        mat_ = std::vector< std::vector<T> >(d - 1);
        for (size_t i = 0; i < mat_.size(); i++)
            mat_[i] = std::vector<T>(mat_.size() - i, 0);
    }
    
    T& operator()(size_t tree, size_t edge)
    {
        return mat_[edge][tree];
    }

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
    
    size_t dim() const 
    {
        return d_;
    }

private:
    size_t d_;
    std::vector< std::vector<T> > mat_;
};



class RVineStructure {
public:
    RVineStructure(const Eigen::MatrixXi& mat)
    {
        d_ = mat.cols();
        trunc_lvl_ = find_trunc_lvl(mat);
        order_ = get_order(mat);
        mat_ = create_struct_mat(mat);
        max_mat_ = compute_max_mat(mat_);
        needed_hfunc1_ = compute_needed_hfunc1(mat_, max_mat_);
        needed_hfunc2_ = compute_needed_hfunc2(mat_, max_mat_);
    }

    RVineMatrix2<size_t> matrix() const {return mat_;}
    RVineMatrix2<size_t> max_matrix() const {return max_mat_;}
    RVineMatrix2<size_t> needed_hfunc1() const {return needed_hfunc1_;}
    RVineMatrix2<size_t> needed_hfunc2() const {return needed_hfunc2_;}

protected:

    size_t find_trunc_lvl(const Eigen::MatrixXi& mat)
    {
        size_t trunc_lvl; ;
        for (trunc_lvl = mat.cols() - 1; trunc_lvl > 0; trunc_lvl--) {
            if (mat(trunc_lvl - 1, 0) != 0)
                break;
        }

        return trunc_lvl;
    }

    std::vector<size_t> get_order(const Eigen::MatrixXi& mat) const
    {
        size_t d = mat.cols();
        std::vector<size_t> order(d);
        for (size_t i = 0; i < d; i++) {
            order[i] = mat(i, d - i - 1);
            if (order[i] != i + 1)
                throw std::runtime_error("mat must be in natural order");
        }
        
        return order;
    }
    
    RVineMatrix2<size_t> create_struct_mat(const Eigen::MatrixXi& mat) const
    {
        RVineMatrix2<size_t> str_mat(d_);
        for (size_t i = 0; i < d_ - 1; i++) {
            for (size_t j = 0; j < d_ - 1 - i; j++)
                str_mat(i, j) = mat(i, j);
        }
        
        return str_mat;
    }
    
    RVineMatrix2<size_t> compute_max_mat(RVineMatrix2<size_t>& mat) const
    {
        RVineMatrix2<size_t> max_mat = mat;
        for (size_t j = 0; j < d_ - 1; j++) {
            for (size_t i = 1; i < d_ - 1 - j; i++) {
                    max_mat(i, j) = std::max(mat(i, j), max_mat(i - 1, j));
            }
        }
    
        return max_mat;
    }
    
    RVineMatrix2<size_t> compute_needed_hfunc1(RVineMatrix2<size_t>& mat,
                                               RVineMatrix2<size_t>& max_mat)
                                               const
    {
        size_t d = mat.dim();
        RVineMatrix2<size_t> needed_hfunc1(d);
 
        for (size_t i = 0; i < d - 2; i++) {
            for (size_t j = 0; j < d - 2 - i; j++) {
                if (mat(i + 1, j) != max_mat(i + 1, j))
                    needed_hfunc1(i, d - max_mat(i + 1, j)) = 1;
            }
        }
        
        return needed_hfunc1;
    }
    
    RVineMatrix2<size_t> compute_needed_hfunc2(RVineMatrix2<size_t>& mat,
                                               RVineMatrix2<size_t>& max_mat)
                                               const
    {
        size_t d = mat.dim();
        RVineMatrix2<size_t> needed_hfunc2(d);
 
        for (size_t j = 0; j < d - 2; j++) {
            for (size_t i = 0; i < d - 2 - j; i++) {
                needed_hfunc2(i, j) = 1;
            }
            if (mat(d - 2 - j, j) == max_mat(d - 2 - j, j))
                needed_hfunc2(d - 2 - j, j) = 1;
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
