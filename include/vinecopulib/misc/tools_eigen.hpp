// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vector>
#include <Eigen/Dense>
#include <boost/property_tree/ptree.hpp>

namespace vinecopulib {

//! Tools for working with Eigen types
namespace tools_eigen
{
    //! An `Eigen::Matrix` containing `bool`s (similar to `Eigen::MatrixXd`).
    typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

    Eigen::Matrix<double, Eigen::Dynamic, 2> swap_cols(
         Eigen::Matrix<double, Eigen::Dynamic, 2> u);
    
    Eigen::VectorXd invert_f(
            const Eigen::VectorXd &x, std::function<Eigen::VectorXd(const Eigen::VectorXd&)> f,
            const double lb = 1e-20,
            const double ub = 1-1e-20,
            int n_iter = 35
    );

    Eigen::Matrix<double, Eigen::Dynamic, 2> expand_grid(
            const Eigen::VectorXd& grid_points
    );

    //! conversion from Eigen::Matrix to boost::property_tree::ptree
    //!
    //! @param matrix the Eigen::Matrix to convert.
    //! @return the corresponding boost::property_tree::ptree.
    template <class T> inline boost::property_tree::ptree matrix_to_ptree(
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix)
    {
        size_t rows = matrix.rows();
        size_t cols = matrix.cols();

        boost::property_tree::ptree output;
        for (size_t i = 0; i < cols; i++) {
            boost::property_tree::ptree col;
            for (size_t j = 0; j < rows; j++) {
                boost::property_tree::ptree cell;
                cell.put_value(matrix(j,i));
                col.push_back(std::make_pair("", cell));
            }
            output.push_back(std::make_pair("", col));
        }

        return output;
    };

    //! conversion from boost::property_tree::ptree to Eigen::Matrix
    //!
    //! @param iroot the boost::property_tree::ptree to convert.
    //! @return the corresponding Eigen::Matrix.
    template <typename T> inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ptree_to_matrix(
            const boost::property_tree::ptree input)
    {

        std::vector<double> vec;
        size_t rows = 0;
        size_t cols = 0;
        for (boost::property_tree::ptree::value_type col : input)
        {
            size_t rows_temp = 0;
            for (boost::property_tree::ptree::value_type cell : col.second)
            {
                rows_temp++;
                vec.push_back(cell.second.get_value<double>());
            }
            if (cols == 0) {
                rows = rows_temp;
                cols++;
            } else if (rows_temp != rows) {
                std::stringstream message;
                message << "column 0 to " << cols-1 << " have " <<
                        rows << " rows, but column" << cols << " has " <<
                        rows_temp << "rows" << std::endl;
                throw std::runtime_error(message.str().c_str());
            } else {
                cols++;
            }
        }

        Eigen::MatrixXd mat = Eigen::MatrixXd::Map(&vec[0], rows, cols);
        return mat.cast <T>();
    };

    Eigen::MatrixXd read_matxd(const char *filename, int max_buffer_size = (int) 1e6);
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> read_matxs(
        const char *filename, int max_buffer_size = (int) 1e6);
}

}
