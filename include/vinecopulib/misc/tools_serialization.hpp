// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <vinecopulib/misc/nlohmann_json.hpp>
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

namespace tools_serialization {

//! conversion from Eigen::Matrix to nlohmann::json
//!
//! @param matrix The Eigen::Matrix to convert.
//! @return the corresponding nlohmann::json.
template<class T>
inline nlohmann::json
matrix_to_json(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix)
{

  nlohmann::json output;
  output["shape"] = { matrix.rows(), matrix.cols() };
  nlohmann::json json_data;
  for (long col = 0; col < matrix.cols(); col++) {
    for (long row = 0; row < matrix.rows(); row++) {
      json_data.push_back(matrix(row, col));
    }
  }
  output["data"] = json_data;

  return output;
}

//! conversion from std::vector to nlohmann::json
//!
//! @param vec The std::vector to convert.
//! @return the corresponding nlohmann::json.
template<class T>
inline nlohmann::json
vector_to_json(const std::vector<T>& vec)
{
  nlohmann::json output = vec;
  return output;
}

//! conversion from nlohmann::json to Eigen::Matrix
//!
//! @param input The nlohmann::json to convert.
//! @return the corresponding Eigen::Matrix.
template<typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
json_to_matrix(const nlohmann::json& input)
{

  size_t rows = input["shape"][0];
  size_t cols = input["shape"][1];

  Eigen::MatrixXd matrix;
  if (!input["data"].is_null()) {
    std::vector<double> vec = input["data"];
    matrix = Eigen::MatrixXd::Map(&vec[0], rows, cols);
  } else {
    matrix = Eigen::MatrixXd();
  }
  return matrix.cast<T>();
}

//! conversion from nlohmann::json to std::vector
//!
//! @param input The nlohmann::json to convert.
//! @return the corresponding std::vector.
template<typename T>
inline std::vector<T>
json_to_vector(const nlohmann::json& input)
{

  std::vector<T> res = input;
  return res;
}

inline nlohmann::json
file_to_json(const std::string& filename)
{
  nlohmann::json output;
  std::ifstream file(filename);
  file >> output;
  return output;
}

inline void
json_to_file(const std::string& filename, const nlohmann::json& json)
{
  std::ofstream file(filename);
  file << json << std::endl;
}
}
}
