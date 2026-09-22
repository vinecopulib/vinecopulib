// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
#include <unordered_map>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_var_types.hpp>

namespace vinecopulib {

namespace {

// The two directions are kept as separate maps rather than a bidirectional
// container: the table is fixed and tiny, and both lookups are on hot paths
// (family names appear in every serialized model).
inline const std::vector<std::pair<BicopFamily, std::string>>&
family_name_table()
{
  static const std::vector<std::pair<BicopFamily, std::string>> table = {
    { BicopFamily::indep, "Independence" },
    { BicopFamily::gaussian, "Gaussian" },
    { BicopFamily::student, "Student" },
    { BicopFamily::clayton, "Clayton" },
    { BicopFamily::gumbel, "Gumbel" },
    { BicopFamily::frank, "Frank" },
    { BicopFamily::joe, "Joe" },
    { BicopFamily::bb1, "BB1" },
    { BicopFamily::bb6, "BB6" },
    { BicopFamily::bb7, "BB7" },
    { BicopFamily::bb8, "BB8" },
    { BicopFamily::tawn, "Tawn" },
    { BicopFamily::tll, "TLL" },
    { BicopFamily::cardioid, "Cardioid" },
    { BicopFamily::wrapped_cauchy, "Wrapped Cauchy" },
    { BicopFamily::von_mises, "von Mises" },
    { BicopFamily::cubic_sections, "Cubic sections" }
  };
  return table;
}

struct BicopFamilyHash
{
  size_t operator()(BicopFamily family) const noexcept
  {
    return std::hash<int>()(static_cast<int>(family));
  }
};

inline const std::unordered_map<BicopFamily, std::string, BicopFamilyHash>&
enum_to_name()
{
  static const std::unordered_map<BicopFamily, std::string, BicopFamilyHash>
    map(family_name_table().begin(), family_name_table().end());
  return map;
}

inline const std::unordered_map<std::string, BicopFamily>&
name_to_enum()
{
  static const std::unordered_map<std::string, BicopFamily> map = [] {
    std::unordered_map<std::string, BicopFamily> m;
    for (const auto& e : family_name_table()) {
      m.emplace(e.second, e.first);
    }
    return m;
  }();
  return map;
}
}

//! @brief Converts a BicopFamily into a string with its name.
//! @param family The family.
inline std::string
get_family_name(BicopFamily family)
{
  return enum_to_name().at(family);
}

//! @brief Converts a string name into a BicopFamily.
//! @param family The family name.
inline BicopFamily
get_family_enum(const std::string& family)
{
  return name_to_enum().at(family);
}

//! @brief Whether a family can model a pair with the given variable types.
//!
//! @details The independence copula and the nonparametric `tll` estimator
//! accept every pair. The other families split by geometry: the linear
//! families (`gaussian`, ..., `tawn`) need two non-circular variables; the
//! binding-density circulas need at least one circular variable; the
//! cylindrical sections copulas need exactly one circular and one linear
//! variable. A circular variable paired with a discrete one is never
//! accepted.
//!
//! @param family The family.
//! @param var_types Two variable types, each `"c"`, `"d"`, or `"a"`.
//! @return Whether `family` is eligible for a pair of these types.
inline bool
family_accepts_var_types(BicopFamily family,
                         const std::vector<std::string>& var_types)
{
  using namespace tools_var_types;
  if (var_types.size() != 2) {
    return false;
  }
  const size_t n_circular =
    is_circular(var_types[0]) + is_circular(var_types[1]);
  if (family == BicopFamily::indep) {
    return true;
  }
  if (family == BicopFamily::tll) {
    // the nonparametric estimator recovers a latent sample of a discrete
    // variable on the normal scale, which has no circular counterpart yet
    return n_circular == 0 || count_discrete(var_types) == 0;
  }
  if (tools_stl::is_member(family, bicop_families::cylindrical)) {
    return n_circular == 1;
  }
  if (tools_stl::is_member(family, bicop_families::circular)) {
    return n_circular > 0;
  }
  return n_circular == 0;
}

//! @brief Filters a family set down to the families eligible for a pair.
//! @param families The candidate families.
//! @param var_types Two variable types, each `"c"`, `"d"`, or `"a"`.
//! @return The members of `families` for which `family_accepts_var_types()`
//!   holds, in their original order.
inline std::vector<BicopFamily>
eligible_families(const std::vector<BicopFamily>& families,
                  const std::vector<std::string>& var_types)
{
  std::vector<BicopFamily> out;
  for (auto fam : families) {
    if (family_accepts_var_types(fam, var_types)) {
      out.push_back(fam);
    }
  }
  return out;
}
}
