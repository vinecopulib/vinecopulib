// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <vinecopulib/misc/tools_optional.hpp>
#include <vinecopulib/misc/tools_stl.hpp>


namespace vinecopulib {

namespace tools_print {

// Printer tags
struct PrintSkip {};        // do nothing
struct PrintDefault {};     // use operator<< on the resolved value
struct PrintYesNo {};       // print "yes"/"no" for bools
struct PrintFamilies {};    // vector<BicopFamily> -> comma-separated names
struct PrintThreads {};     // size_t -> max(1, value)
struct PrintWeights {};     // Eigen::VectorXd -> "yes"/"no" by size>0

// Generic
template <class Os, class Opt>
void print_field(Os& os, const char* label, const Opt& opt, PrintDefault) {
  os << label << opt.value() << '\n';
}

// Intentionally empty
template <class Os, class Opt>
void print_field(Os&, const char*, const Opt&, PrintSkip) {

}

// Booleans as yes/no
template <class Os, class Opt>
void print_field(Os& os, const char* label, const Opt& opt, PrintYesNo) {
  os << label << (opt.value() ? "yes" : "no") << '\n';
}

// Families: vector<BicopFamily>
template <class Os>
void print_field(Os& os, const char* label,
                 const optional::optional<std::vector<BicopFamily>>& opt,
                 PrintFamilies) {
  os << label;
  const auto& fs = opt.value();
  for (size_t j = 0; j < fs.size(); ++j) {
    if (j) os << ", ";
    os << get_family_name(fs[j]);
  }
  os << '\n';
}

// Threads
template <class Os, class Opt>
void print_field(Os& os, const char* label, const Opt& opt, PrintThreads) {
  auto n = opt.value();
  os << label << (n == 0 ? size_t{1} : n) << '\n';
}

// Weights: Eigen::VectorXd -> "yes"/"no" by size>0
template <class Os>
void print_field(Os& os, const char* label,
                 const optional::optional<Eigen::VectorXd>& opt,
                 PrintWeights) {
  const auto& w = opt.value();
  os << label << (w.size() > 0 ? "yes" : "no") << '\n';
}


// Function to format vectors of strings like a DataFrame and return a
// stringstream
inline std::stringstream
dataframe_to_string(const std::vector<std::vector<std::string>>& columns)
{

  // TODO: Check if all vectors have the same length

  // Determine maximum column width for pretty printing
  std::vector<size_t> col_widths;
  for (const auto& col : columns) {
    size_t max_width = 0;
    for (const auto& item : col) {
      max_width = std::max(max_width, item.size());
    }
    col_widths.push_back(max_width);
  }

  std::stringstream ss;
  // Assuming all columns have the same number of rows
  size_t num_rows = columns[0].size();
  for (size_t row = 0; row < num_rows; ++row) {
    for (size_t col = 0; col < columns.size(); ++col) {
      ss << std::setw(static_cast<int>(col_widths[col])) << columns[col][row]
         << " ";
    }
    ss << std::endl;
  }

  return ss;
}

} // namespace tools_print
} // namespace vinecopulib
