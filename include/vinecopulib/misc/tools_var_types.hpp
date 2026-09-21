// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <string>
#include <vector>

namespace vinecopulib {

//! @brief Predicates on the `var_types` tokens.
//!
//! A variable type is one of the strings `"c"` (continuous), `"d"`
//! (discrete), or `"a"` (continuous circular). The tokens are compared only
//! here, so that code paths ask the question they mean: whether a variable
//! is discrete, whether it is circular, or whether it is continuous at all.
namespace tools_var_types {

//! the continuous linear type.
inline const std::string&
continuous()
{
  static const std::string type = "c";
  return type;
}

//! the discrete type.
inline const std::string&
discrete()
{
  static const std::string type = "d";
  return type;
}

//! the continuous circular (angular) type.
inline const std::string&
circular()
{
  static const std::string type = "a";
  return type;
}

//! whether `type` is the discrete type.
inline bool
is_discrete(const std::string& type)
{
  return type == discrete();
}

//! whether `type` is the circular type.
inline bool
is_circular(const std::string& type)
{
  return type == circular();
}

//! whether `type` is a continuous type, linear or circular.
inline bool
is_continuous(const std::string& type)
{
  return !is_discrete(type);
}

//! whether `type` is the continuous linear type.
inline bool
is_linear(const std::string& type)
{
  return type == continuous();
}

//! whether `type` is one of the three known types.
inline bool
is_valid(const std::string& type)
{
  return is_linear(type) || is_discrete(type) || is_circular(type);
}

//! whether no entry of `types` is discrete.
inline bool
all_continuous(const std::vector<std::string>& types)
{
  for (const auto& t : types) {
    if (is_discrete(t)) {
      return false;
    }
  }
  return true;
}

//! whether every entry of `types` is the continuous linear type.
inline bool
all_linear(const std::vector<std::string>& types)
{
  for (const auto& t : types) {
    if (!is_linear(t)) {
      return false;
    }
  }
  return true;
}

//! whether every entry of `types` is discrete.
inline bool
all_discrete(const std::vector<std::string>& types)
{
  for (const auto& t : types) {
    if (!is_discrete(t)) {
      return false;
    }
  }
  return true;
}

//! whether any entry of `types` is circular.
inline bool
any_circular(const std::vector<std::string>& types)
{
  for (const auto& t : types) {
    if (is_circular(t)) {
      return true;
    }
  }
  return false;
}

//! the number of discrete entries of `types`.
inline size_t
count_discrete(const std::vector<std::string>& types)
{
  size_t n = 0;
  for (const auto& t : types) {
    n += is_discrete(t);
  }
  return n;
}

//! `types` with every discrete entry replaced by the continuous linear type;
//! circular entries are kept.
inline std::vector<std::string>
as_continuous(std::vector<std::string> types)
{
  for (auto& t : types) {
    if (is_discrete(t)) {
      t = continuous();
    }
  }
  return types;
}

//! a vector of `n` continuous linear types, the default `var_types`.
inline std::vector<std::string>
all_continuous_types(size_t n)
{
  return std::vector<std::string>(n, continuous());
}

} // namespace tools_var_types
} // namespace vinecopulib
