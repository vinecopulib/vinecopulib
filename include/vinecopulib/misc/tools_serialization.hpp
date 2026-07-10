// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <stdexcept>
#include <vector>

// GCC 11 diagnoses a fallthrough in the vendored CBOR parser when its stream
// overload is instantiated.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
#include <vinecopulib/misc/nlohmann_json.hpp>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

//! @brief Encoding used for model persistence.
//!
//! JSON is the human-readable interchange format. CBOR is a binary encoding
//! of the same logical `nlohmann::json` representation.
enum class SerializationFormat
{
  json,
  cbor
};

namespace tools_serialization {

inline const char*
serialization_format_name(const SerializationFormat format)
{
  switch (format) {
    case SerializationFormat::json:
      return "JSON";
    case SerializationFormat::cbor:
      return "CBOR";
    default:
      throw std::invalid_argument("unsupported serialization format");
  }
}

inline SerializationFormat
serialization_format_from_filename(const std::string& filename)
{
  const std::string extension = ".cbor";
  if (filename.size() >= extension.size() &&
      filename.compare(
        filename.size() - extension.size(), extension.size(), extension) == 0) {
    return SerializationFormat::cbor;
  }
  return SerializationFormat::json;
}

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

//! conversion from vinecopulib::TriangularArray to nlohmann::json
//!
//! @param array The vinecopulib::TriangularArray to convert.
//! @return the corresponding nlohmann::json.
template<class T>
inline nlohmann::json
triangular_array_to_json(const TriangularArray<T>& array)
{
  nlohmann::json output;
  size_t d = array.get_dim();
  size_t trunc_lvl = array.get_trunc_lvl();
  output["d"] = d;
  output["t"] = trunc_lvl;

  nlohmann::json json_data;
  for (size_t i = 0; i < std::min(d - -1, trunc_lvl); i++) {
    nlohmann::json row;
    for (size_t j = 0; j < d - 1 - i; j++) {
      row.push_back(array(i, j));
    }
    json_data.push_back(row);
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

//! conversion from nlohmann::json to vinecopulib::TriangularArray
//!
//! @param input The nlohmann::json to convert.
//! @return the corresponding vinecopulib::TriangularArray
template<typename T>
inline TriangularArray<T>
json_to_triangular_array(const nlohmann::json& input)
{

  std::vector<std::vector<T>> vec = input["data"];
  return TriangularArray<T>(vec);
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
file_to_json(const std::string& filename, const SerializationFormat format)
{
  const auto format_name = serialization_format_name(format);
  auto mode = std::ios::in;
  if (format == SerializationFormat::cbor) {
    mode |= std::ios::binary;
  }
  std::ifstream file(filename, mode);
  if (!file.is_open()) {
    throw std::runtime_error("could not open " + filename + " for reading");
  }

  try {
    if (format == SerializationFormat::cbor) {
      return nlohmann::json::from_cbor(file);
    }
    nlohmann::json output;
    file >> output;
    return output;
  } catch (const nlohmann::json::exception& exception) {
    throw std::runtime_error("failed to parse " + std::string(format_name) +
                             " file " + filename + ": " + exception.what());
  }
}

inline void
json_to_file(const std::string& filename,
             const nlohmann::json& json,
             const SerializationFormat format)
{
  serialization_format_name(format);
  auto mode = std::ios::out;
  if (format == SerializationFormat::cbor) {
    mode |= std::ios::binary;
  }
  std::ofstream file(filename, mode);
  if (!file.is_open()) {
    throw std::runtime_error("could not open " + filename + " for writing");
  }

  if (format == SerializationFormat::cbor) {
    nlohmann::json::to_cbor(json, file);
  } else {
    file << json << std::endl;
  }
  file.flush();
  if (!file) {
    throw std::runtime_error("failed to write " + filename);
  }
}

inline nlohmann::json
file_to_json(const std::string& filename)
{
  return file_to_json(filename, serialization_format_from_filename(filename));
}

inline void
json_to_file(const std::string& filename, const nlohmann::json& json)
{
  json_to_file(filename, json, serialization_format_from_filename(filename));
}
}
}
