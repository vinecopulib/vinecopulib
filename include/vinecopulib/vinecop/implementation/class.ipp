// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

#include <stdexcept>

namespace vinecopulib {

//! @brief creates a D-vine on `d` variables with all pair-copulas set to
//! independence.
//! @param d the dimension (= number of variables) of the model.
inline Vinecop::Vinecop(const size_t d)
  : Vinecop(RVineStructure(tools_stl::seq_int(1, d), static_cast<size_t>(0)))
{}

//! @brief creates a vine copula with structure specified by an RVineStructure
//! object; all pair-copulas are set to independence.
//! @param structure an RVineStructure object representing the structure of
//! the vine.
inline Vinecop::Vinecop(const RVineStructure& structure)
  : d_(structure.get_dim())
  , vine_struct_(structure)
  // pair_copulas_ empty = everything independence
  , threshold_(0.0)
  , loglik_(NAN)
{
  set_continuous_var_types();
}

//! @brief creates a vine copula with structure specified by an R-vine matrix;
//! all pair-copulas are set to independence.
//! @param matrix an R-vine matrix.
//! @param check whether to check if `matrix` is a valid R-vine
//!     matrix.
inline Vinecop::Vinecop(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
  const bool check)
  : Vinecop(RVineStructure(matrix, check))
{}

//! @brief creates a vine copula with structure specified by an R-vine matrix;
//! all pair-copulas are set to independence.
//! @param order the order of the variables in the vine structure, see
//! RVineStructure's corresponding constructor.
//! @param struct_array a triangular array object specifying the vine structure,
//! see RVineStructure's corresponding constructor.
//! @param check whether `order` and `struct_array` shall be checked
//! for validity.
inline Vinecop::Vinecop(const std::vector<size_t>& order,
                        const TriangularArray<size_t>& struct_array,
                        const bool check)
  : Vinecop(RVineStructure(order, struct_array, false, check))
{}

//! @brief creates an arbitrary vine copula model.
//! @param pair_copulas Bicop objects specifying the pair-copulas, see
//!     make_pair_copula_store().
//! @param structure an RVineStructure object specifying the vine structure.
inline Vinecop::Vinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
                        const RVineStructure& structure)
  : Vinecop(structure)
{
  check_pair_copulas_rvine_structure(pair_copulas);
  pair_copulas_ = pair_copulas;
  vine_struct_.truncate(pair_copulas.size());
}

//! @brief creates an arbitrary vine copula model.
//! @param pair_copulas Bicop objects specifying the pair-copulas, see
//!     make_pair_copula_store().
//! @param matrix an R-vine matrix specifying the vine structure.
//! @param check whether to check if `matrix` is a valid R-vine
//!     matrix.
inline Vinecop::Vinecop(
  const std::vector<std::vector<Bicop>>& pair_copulas,
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
  const bool check)
  : Vinecop(pair_copulas, RVineStructure(matrix, check))
{}

//! @brief creates an arbitrary vine copula model.
//! @param pair_copulas Bicop objects specifying the pair-copulas, see
//!     make_pair_copula_store().
//! @param order the order of the variables in the vine structure, see
//! RVineStructure's corresponding constructor.
//! @param struct_array a triangular array object specifying the vine structure,
//! see the corresponding constructor of RVineStructure.
//! @param check whether `order` and `struct_array` shall be checked
//! for validity.
inline Vinecop::Vinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
                        const std::vector<size_t>& order,
                        const TriangularArray<size_t>& struct_array,
                        const bool check)
  : Vinecop(pair_copulas, RVineStructure(order, struct_array, false, check))
{}

//! @brief creates from a boost::property_tree::ptree object
//! @param input the boost::property_tree::ptree object to convert from
//! (see to_ptree() for the structure of the input).
//! @param check whether to check if the `"structure"` node represents
//!      a valid R-vine structure.
inline Vinecop::Vinecop(const boost::property_tree::ptree input,
                        const bool check)
{
  vine_struct_ = RVineStructure(input.get_child("structure"), check);
  d_ = static_cast<size_t>(vine_struct_.get_dim());

  boost::property_tree::ptree pcs_node = input.get_child("pair copulas");
  for (size_t tree = 0; tree < d_ - 1; ++tree) {
    boost::property_tree::ptree tree_node;
    try {
      tree_node = pcs_node.get_child("tree" + std::to_string(tree));
    } catch (...) {
      break; // vine was truncated, no more trees to parse
    }
    // reserve space for pair copulas of this tree
    pair_copulas_.resize(tree + 1);
    pair_copulas_[tree].resize(d_ - tree - 1);

    for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
      boost::property_tree::ptree pc_node =
        tree_node.get_child("pc" + std::to_string(edge));
      pair_copulas_[tree][edge] = Bicop(pc_node);
    }
  }

  // try block for backwards compatibility
  try {
    var_types_ = tools_serialization::ptree_to_vector<std::string>(
      input.get_child("var_types"));
    nobs_ = input.get<size_t>("nobs_");
    threshold_ = input.get<double>("threshold");
    loglik_ = input.get<double>("loglik");
  } catch (...) {
  }
}

//! @brief creates from a JSON file.
//! @param filename the name of the JSON file to read (see to_ptree() for the
//! structure of the file).
//! @param check whether to check if the `"structure"` node of the input
//! represents
//!      a valid R-vine structure.
inline Vinecop::Vinecop(const std::string filename, const bool check)
  : Vinecop(tools_serialization::json_to_ptree(filename.c_str()), check)
{}

//! @brief constructs a vine copula model from data by creating a model and
//! calling select().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param structure an RVineStructure object specifying the vine structure.
//! @param controls see FitControlsVinecop.
inline Vinecop::Vinecop(const Eigen::MatrixXd& data,
                        const RVineStructure& structure,
                        FitControlsVinecop controls)
  : Vinecop(structure)
{
  nobs_ = data.rows();
  check_enough_data(data);
  if (static_cast<size_t>(data.cols()) != d_) {
    throw std::runtime_error("data and structure have "
                             "incompatible dimensions.");
  }
  check_weights_size(controls.get_weights(), data);
  select(data, controls);
}

//! @brief constructs a vine copula model from data by creating a model and
//! calling select().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param matrix either an empty matrix (default) or an R-vine structure
//!     matrix, see select().
//! @param controls see FitControlsVinecop.
//! @param check whether to check if `matrix` is a valid R-vine
//!     matrix.
inline Vinecop::Vinecop(
  const Eigen::MatrixXd& data,
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
  FitControlsVinecop controls,
  const bool check)
  : Vinecop(data, RVineStructure(matrix, check), controls)
{}

//! @brief constructs a vine copula model from data by creating a model and
//! calling select().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param order the order of the variables in the vine structure, see
//! the corresponding constructor of RVineStructure.
//! @param struct_array a triangular array object specifying the vine structure,
//! see the corresponding constructor of RVineStructure.
//! @param controls see FitControlsVinecop.
//! @param check whether `order` and `struct_array` shall be checked
//! for validity.
inline Vinecop::Vinecop(const Eigen::MatrixXd& data,
                        const std::vector<size_t>& order,
                        const TriangularArray<size_t>& struct_array,
                        FitControlsVinecop controls,
                        const bool check)
  : Vinecop(data, RVineStructure(order, struct_array, false, check), controls)
{}

//! @brief constructs a vine copula model from data by creating a model and
//! calling select().
//!
//! @param data an \f$ n \times d \f$ matrix of observations.
//! @param controls see FitControlsVinecop.
inline Vinecop::Vinecop(const Eigen::MatrixXd& data,
                        const FitControlsVinecop& controls)
  : Vinecop(data.cols())
{
  nobs_ = data.rows();
  check_enough_data(data);
  check_weights_size(controls.get_weights(), data);
  select(data, controls);
}

//! @brief converts the copula into a boost::property_tree::ptree object.
//!
//! The `ptree` object contains two nodes : `"structure"` for the vine
//! structure, which itself contains nodes `"array"` for the structure
//! triangular array and `"order"` for the order vector, and `"pair copulas"`.
//! The former two encode the R-Vine structure and the latter is a list of
//! child nodes for the trees (`"tree1"`, `"tree2"`, etc), each containing
//! a list of child nodes for the edges (`"pc1"`, `"pc2"`, etc).
//! See Bicop::to_ptree() for the encoding of pair-copulas.
//!
//! @return the boost::property_tree::ptree object containing the copula.
inline boost::property_tree::ptree
Vinecop::to_ptree() const
{
  boost::property_tree::ptree pair_copulas;
  for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
    boost::property_tree::ptree tree_node;
    for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
      tree_node.add_child("pc" + std::to_string(edge),
                          pair_copulas_[tree][edge].to_ptree());
    }
    pair_copulas.add_child("tree" + std::to_string(tree), tree_node);
  }

  boost::property_tree::ptree output;
  output.add_child("pair copulas", pair_copulas);
  auto structure_node = vine_struct_.to_ptree();
  output.add_child("structure", structure_node);
  output.add_child("var_types",
                   tools_serialization::vector_to_ptree(var_types_));
  output.put("nobs_", nobs_);
  output.put("threshold", threshold_);
  output.put("loglik", loglik_);

  return output;
}

//! @brief write the copula object into a JSON file.
//!
//! See to_ptree() for the structure of the file.
//! @param filename the name of the file to write.
inline void
Vinecop::to_json(const std::string filename) const
{
  boost::property_tree::write_json(filename.c_str(), this->to_ptree());
}

//! @brief initializes object for storing pair copulas.
//!
//! @param d dimension of the vine copula.
//! @param trunc_lvl a truncation level (optional).
//! @return A nested vector such that `pc_store[t][e]` contains a Bicop.
//!     object for the pair copula corresponding to tree `t` and edge `e`.
inline std::vector<std::vector<Bicop>>
Vinecop::make_pair_copula_store(const size_t d, const size_t trunc_lvl)
{
  return tools_select::VinecopSelector::make_pair_copula_store(d, trunc_lvl);
}

//! @brief automatically fits and selects a vine copula model.
//!
//! @details `select()` behaves differently depending on its current truncation
//! level and the truncation level specified in the controls, respectively
//! called `trunc_lvl` and `controls.trunc_lvl` in what follows.
//! Essentially, `controls.trunc_lvl` defines the object's truncation level
//! after calling `select()`:
//!   - If `controls.trunc_lvl <= trunc_lvl`, the families and parameters for
//!     all pairs in trees smaller or equal to `controls.trunc_lvl`
//!     are selected, using the current structure.
//!   - If `controls.trunc_lvl > trunc_lvl`, `select()` behaves as above for
//!     all trees that are smaller or equal to `trunc_lvl`, and then it selects
//!     the structure for higher trees along with the families and parameters.
//!     This includes the case where `trunc_lvl = 0`, namely where the
//!     structure is fully unspecified.
//!
//! Selection of the structure is performed using the algorithm of
//! Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
//! *Selecting and estimating regular vine copulae and application to
//! financial returns.* Computational Statistics & Data Analysis, 59 (1),
//! 52-69.
//!
//! When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times d \f$ block contains realizations of
//! \f$ F_Y(Y), F_X(X) \f$; the second \f$ n \times d \f$ block contains
//! realizations of \f$ F_Y(Y^-), F_X(X^-), ... \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_Y(Y^-) = F_Y(Y - 1) \f$. Continuous variables in the second block can
//! be omitted.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls the controls to the algorithm (see FitControlsVinecop).
inline void
Vinecop::select(Eigen::MatrixXd data, const FitControlsVinecop& controls)
{
  check_data(data);
  data = collapse_data(data);

  tools_select::VinecopSelector selector(
    data, vine_struct_, controls, var_types_);
  if (controls.needs_sparse_select()) {
    selector.sparse_select_all_trees(data);
  } else {
    selector.select_all_trees(data);
  }
  finalize_fit(selector);
}

//! @brief automatically fits and selects a vine copula model.
//!
//! Selection of the structure is performed using the algorithm of
//! Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
//! *Selecting and estimating regular vine copulae and application to
//! financial returns.* Computational Statistics & Data Analysis, 59 (1),
//! 52-69.
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times d \f$ block contains realizations of
//! \f$ F_Y(Y), F_X(X) \f$; the second \f$ n \times d \f$ block contains
//! realizations of \f$ F_Y(Y^-), F_X(X^-), ... \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_Y(Y^-) = F_Y(Y - 1) \f$. Continuous variables in the second block can
//! be omitted.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls the controls to the algorithm (see FitControlsVinecop).
inline void
Vinecop::select_all(Eigen::MatrixXd data, const FitControlsVinecop& controls)
{
  vine_struct_ =
    RVineStructure(tools_stl::seq_int(1, d_), static_cast<size_t>(0));
  select(data, controls);
}

//! @brief automatically selects all pair-copula families and fits all
//! parameters.
//!
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times d \f$ block contains realizations of
//! \f$ F_Y(Y), F_X(X) \f$; the second \f$ n \times d \f$ block contains
//! realizations of \f$ F_Y(Y^-), F_X(X^-), ... \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_Y(Y^-) = F_Y(Y - 1) \f$. Continuous variables in the second block can
//! be omitted.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls the controls to the algorithm (see FitControlsVinecop).
inline void
Vinecop::select_families(Eigen::MatrixXd data,
                         const FitControlsVinecop& controls)
{
  controls.set_trunc_lvl(vine_struct_.get_trunc_lvl());
  select(data, controls);
}

//! @name Getters
//! @{

//! @brief extracts a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline Bicop
Vinecop::get_pair_copula(const size_t tree, const size_t edge) const
{
  if (tree > d_ - 2) {
    std::stringstream message;
    message << "tree index out of bounds" << std::endl
            << "allowed: 0, ..., " << d_ - 2 << std::endl
            << "actual: " << tree << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  if (edge > d_ - tree - 2) {
    std::stringstream message;
    message << "edge index out of bounds" << std::endl
            << "allowed: 0, ..., " << d_ - tree - 2 << std::endl
            << "actual: " << edge << std::endl
            << "tree level: " << tree << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  if (tree >= pair_copulas_.size()) {
    // vine is truncated
    return Bicop();
  }
  return pair_copulas_[tree][edge];
}

//! @brief extracts all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<Bicop>>
Vinecop::get_all_pair_copulas() const
{
  return pair_copulas_;
}

//! @brief extracts the family of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline BicopFamily
Vinecop::get_family(const size_t tree, const size_t edge) const
{
  return get_pair_copula(tree, edge).get_family();
}

//! @brief extracts the families of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<BicopFamily>>
Vinecop::get_all_families() const
{
  std::vector<std::vector<BicopFamily>> families(pair_copulas_.size());
  for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
    families[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      families[tree][edge] = get_family(tree, edge);
    }
  }

  return families;
}

//! @brief extracts the rotation of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline int
Vinecop::get_rotation(const size_t tree, const size_t edge) const
{
  return get_pair_copula(tree, edge).get_rotation();
}

//! @brief extracts the rotations of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<int>>
Vinecop::get_all_rotations() const
{
  std::vector<std::vector<int>> rotations(pair_copulas_.size());
  for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
    rotations[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      rotations[tree][edge] = get_rotation(tree, edge);
    }
  }

  return rotations;
}

//! @brief extracts the parameters of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline Eigen::MatrixXd
Vinecop::get_parameters(const size_t tree, const size_t edge) const
{
  return get_pair_copula(tree, edge).get_parameters();
}

//! @brief extracts the Kendall's \f$ tau \f$ of a pair copula.
//!
//! @param tree tree index (starting with 0).
//! @param edge edge index (starting with 0).
inline double
Vinecop::get_tau(const size_t tree, const size_t edge) const
{
  return get_pair_copula(tree, edge).get_tau();
}

inline size_t
Vinecop::get_trunc_lvl() const
{
  return vine_struct_.get_trunc_lvl();
}

//! @brief extracts the parameters of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<Eigen::MatrixXd>>
Vinecop::get_all_parameters() const
{
  std::vector<std::vector<Eigen::MatrixXd>> parameters(pair_copulas_.size());
  for (size_t tree = 0; tree < parameters.size(); ++tree) {
    parameters[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      parameters[tree][edge] = get_parameters(tree, edge);
    }
  }

  return parameters;
}

//! @brief extracts the Kendall's \f$ tau \f$s of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<double>>
Vinecop::get_all_taus() const
{
  std::vector<std::vector<double>> taus(pair_copulas_.size());
  for (size_t tree = 0; tree < taus.size(); ++tree) {
    taus[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      taus[tree][edge] = get_tau(tree, edge);
    }
  }

  return taus;
}

//! @brief extracts the dimension of the vine copula model.
inline size_t
Vinecop::get_dim() const
{
  return d_;
}

//! @brief extracts the order vector of the vine copula model.
inline std::vector<size_t>
Vinecop::get_order() const
{
  return vine_struct_.get_order();
}

//! @brief extracts the structure matrix of the vine copula model.
inline RVineStructure
Vinecop::get_rvine_structure() const
{
  return vine_struct_;
}

//! @brief extracts the structure matrix of the vine copula model.
inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
Vinecop::get_matrix() const
{
  return vine_struct_.get_matrix();
}

//! @brief extracts the above diagonal coefficients of the vine copula model.
//! @param natural_order whether indices correspond to natural order.
inline TriangularArray<size_t>
Vinecop::get_struct_array(bool natural_order) const
{
  return vine_struct_.get_struct_array(natural_order);
}

//! @brief extracts the log-likelihood (throws an error if model has not been
//! fitted to data).
inline double
Vinecop::get_loglik() const
{
  check_fitted();
  return loglik_;
}

//! @brief extracts the number of observations used for the fit.
//!
//! The function throws an error if model has not been fitted to data.
inline size_t
Vinecop::get_nobs() const
{
  check_fitted();
  return nobs_;
}

//! @brief extracts the AIC.
//!
//! The function throws an error if model has not been fitted to data.
inline double
Vinecop::get_aic() const
{
  check_fitted();
  return -2 * loglik_ + 2 * get_npars();
}

//! @brief extracts the BIC.
//!
//! The function throws an error if model has not been fitted to data.
inline double
Vinecop::get_bic() const
{
  check_fitted();
  return -2 * loglik_ + get_npars() * std::log(nobs_);
}

//! @brief extracts the log-likelihood.
//!
//! The function throws an error if model has not been fitted to data.
inline double
Vinecop::get_mbicv(const double psi0) const
{
  check_fitted();
  return -2 * loglik_ + this->calculate_mbicv_penalty(nobs_, psi0);
}

//! computes the penalty term for mBICV
inline double
Vinecop::calculate_mbicv_penalty(const size_t nobs, const double psi0) const
{
  if (!(psi0 > 0.0) | !(psi0 < 1.0)) {
    throw std::runtime_error("psi0 must be in the interval (0, 1)");
  }
  auto all_fams = get_all_families();
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> non_indeps(d_ - 1);
  non_indeps.setZero();
  for (size_t t = 0; t < d_ - 1; t++) {
    if (t == all_fams.size()) {
      break;
    }
    for (size_t e = 0; e < d_ - 1 - t; e++) {
      if (all_fams[t][e] != BicopFamily::indep) {
        non_indeps(t)++;
      }
    }
  }
  auto sq0 = tools_stl::seq_int(1, d_ - 1);
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> sq(d_ - 1);
  auto psis = Eigen::VectorXd(d_ - 1);
  for (size_t i = 0; i < d_ - 1; i++) {
    sq(i) = sq0[i];
    psis(i) = std::pow(psi0, sq0[i]);
  }
  double npars = this->get_npars();
  double log_prior = (non_indeps.cast<double>().array() * psis.array().log() +
                      (d_ - non_indeps.array() - sq.array()).cast<double>() *
                        (1 - psis.array()).log())
                       .sum();

  return std::log(nobs) * npars - 2 * log_prior;
}

//! @brief extracts the threshold (usually zero except `select_threshold ==
//! TRUE` in `FitControlsVinecop()`).
inline double
Vinecop::get_threshold() const
{
  return threshold_;
}

//! @brief sets variable types.
//! @param var_types a vector specifying the types of the variables,
//!   e.g., `{"c", "d"}` means first varible continuous, second discrete.
inline void
Vinecop::set_var_types(const std::vector<std::string>& var_types)
{
  check_var_types(var_types);
  set_var_types_internal(var_types);
}

inline void
Vinecop::check_var_types(const std::vector<std::string>& var_types) const
{
  std::stringstream msg;
  if (var_types.size() > d_) {
    msg << "more var_types (" << var_types.size() << ")"
        << "than variables (" << d_ << ")-" << std::endl;
    throw std::runtime_error(msg.str());
  }
  for (auto t : var_types) {
    if (!tools_stl::is_member(t, { "c", "d" })) {
      msg << "variable type must be 'c' or 'd' (not '" << t << "')."
          << std::endl;
      throw std::runtime_error(msg.str());
    }
  }
}

//! @brief sets variable types.
//! @param var_types a vector specifying the types of the variables,
//!   e.g., `{"c", "d"}` means first varible continuous, second discrete.
inline void
Vinecop::set_var_types_internal(const std::vector<std::string>& var_types) const
{
  var_types_ = var_types;
  if (pair_copulas_.size() == 0) {
    return;
  }

  // set new var_types for all pair-copulas
  std::vector<std::string> natural_types(d_), pair_types(2);
  for (size_t j = 0; j < d_; ++j) {
    natural_types[j] = var_types[vine_struct_.get_order()[j] - 1];
  }
  // we set the first tree explicitly and deduce later trees
  for (size_t e = 0; e < d_ - 1; ++e) {
    pair_types[0] = natural_types[e];
    pair_types[1] = natural_types[vine_struct_.struct_array(0, e, true) - 1];
    pair_copulas_[0][e].set_var_types(pair_types);
  }

  for (size_t t = 1; t < pair_copulas_.size(); ++t) {
    for (size_t e = 0; e < d_ - t - 1; ++e) {
      size_t m = vine_struct_.min_array(t, e);
      pair_types[0] = pair_copulas_[t - 1][e].get_var_types()[0];
      if (m == vine_struct_.struct_array(t, e, true)) {
        pair_types[1] = pair_copulas_[t - 1][m - 1].get_var_types()[0];
      } else {
        pair_types[1] = pair_copulas_[t - 1][m - 1].get_var_types()[1];
      }
      pair_copulas_[t][e].set_var_types(pair_types);
    }
  }
}

//! @brief extracts the variable types.
std::vector<std::string>
Vinecop::get_var_types() const
{
  return var_types_;
}

//! @}

//! @brief calculates the density function of the vine copula model.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::VectorXd
Vinecop::pdf(Eigen::MatrixXd u, const size_t num_threads) const
{
  check_data(u);
  u = collapse_data(u);

  size_t d = d_;
  size_t n = u.rows();

  // info about the vine structure (reverse rows (!) for more natural indexing)
  size_t trunc_lvl = vine_struct_.get_trunc_lvl();
  std::vector<size_t> order;
  TriangularArray<size_t> natural_array, min_array, needed_hfunc1,
    needed_hfunc2;
  if (trunc_lvl > 0) {
    order = vine_struct_.get_order();
    natural_array = vine_struct_.get_struct_array(true);
    min_array = vine_struct_.get_min_array();
    needed_hfunc1 = vine_struct_.get_needed_hfunc1();
    needed_hfunc2 = vine_struct_.get_needed_hfunc2();
  }
  auto disc_cols = tools_select::get_disc_cols(var_types_);

  // initial value must be 1.0 for multiplication
  Eigen::VectorXd pdf = Eigen::VectorXd::Constant(u.rows(), 1.0);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects (all data must be in (0, 1))
    Eigen::MatrixXd hfunc1(b.size, d);
    hfunc1.setZero();
    Eigen::MatrixXd hfunc2 = hfunc1;
    Eigen::MatrixXd hfunc1_sub = hfunc1;
    Eigen::MatrixXd hfunc2_sub = hfunc1;
    Eigen::MatrixXd u_e(b.size, 4);
    Eigen::MatrixXd u_sub(b.size, 4);

    // fill first row of hfunc2 matrix with evaluation points;
    // points have to be reordered to correspond to natural order

    for (size_t j = 0; j < d; ++j) {
      hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
      if (var_types_[order[j] - 1] == "d") {
        hfunc2_sub.col(j) =
          u.block(b.begin, d_ + disc_cols[order[j] - 1], b.size, 1);
      } else {
        hfunc2_sub.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
      }
    }

    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      tools_interface::check_user_interrupt(n * d > 1e5);
      for (size_t edge = 0; edge < d - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        u_e.col(0) = hfunc2.col(edge);
        u_e.col(2) = hfunc2_sub.col(edge);
        size_t m = min_array(tree, edge);
        if (m == natural_array(tree, edge)) {
          u_e.col(1) = hfunc2.col(m - 1);
          u_e.col(3) = hfunc2_sub.col(m - 1);
        } else {
          u_e.col(1) = hfunc1.col(m - 1);
          u_e.col(3) = hfunc1_sub.col(m - 1);
        }

        Bicop edge_copula = get_pair_copula(tree, edge);
        pdf.segment(b.begin, b.size) =
          pdf.segment(b.begin, b.size).cwiseProduct(edge_copula.pdf(u_e));

        // h-functions are only evaluated if needed in next step
        auto var_types = edge_copula.get_var_types();
        if (needed_hfunc1(tree, edge)) {
          hfunc1.col(edge) = edge_copula.hfunc1(u_e);
          if (var_types[1] == "d") {
            u_sub = u_e;
            u_sub.col(1) = u_sub.col(3);
            hfunc1_sub.col(edge) = edge_copula.hfunc1(u_sub);
          }
        }
        if (needed_hfunc2(tree, edge)) {
          hfunc2.col(edge) = edge_copula.hfunc2(u_e);
          if (var_types[0] == "d") {
            u_sub = u_e;
            u_sub.col(0) = u_sub.col(2);
            hfunc2_sub.col(edge) = edge_copula.hfunc2(u_sub);
          }
        }
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();
  }

  return pdf;
}

//! @brief calculates the cumulative distribution of the vine copula model.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param N integer for the number of quasi-random numbers to draw
//! to evaluate the distribution (default: 1e4).
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will generate `n` samples concurrently in
//!   `num_threads` batches.
//! @param seeds seeds to scramble the quasi-random numbers; if empty (default),
//!   the random number quasi-generator is seeded randomly.
inline Eigen::VectorXd
Vinecop::cdf(const Eigen::MatrixXd& u,
             const size_t N,
             const size_t num_threads,
             std::vector<int> seeds) const
{
  if (d_ > 21201) {
    std::stringstream message;
    message << "cumulative distribution available for models of "
            << "dimension 21201 or less. This model's dimension: " << d_
            << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  check_data(u);

  // Simulate N quasi-random numbers from the vine model
  auto u_sim = simulate(N, true, num_threads, seeds);

  size_t n = u.rows();
  Eigen::VectorXd vine_distribution(n);
  Eigen::ArrayXXd x(N, 1);
  Eigen::RowVectorXd temp(d_);
  for (size_t i = 0; i < n; i++) {
    tools_interface::check_user_interrupt(i % 1000 == 0);
    temp = u.block(i, 0, 1, d_);
    x = (u_sim.rowwise() - temp).rowwise().maxCoeff().array();
    vine_distribution(i) = (x <= 0.0).count();
  }
  return vine_distribution / static_cast<double>(N);
}

//! @brief simulates from a vine copula model, see inverse_rosenblatt().
//!
//! @details Simulated data is always a continous \f$ n \times d \f$ matrix.
//!
//! @param n number of observations.
//! @param qrng set to true for quasi-random numbers.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will generate `n` samples concurrently in
//!   `num_threads` batches.
//! @param seeds seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @return An \f$ n \times d \f$ matrix of samples from the copula model.
inline Eigen::MatrixXd
Vinecop::simulate(const size_t n,
                  const bool qrng,
                  const size_t num_threads,
                  const std::vector<int>& seeds) const
{
  auto u = tools_stats::simulate_uniform(n, d_, qrng, seeds);
  // inverse_rosenblatt() only works for continous models
  auto actual_types = var_types_;
  set_continuous_var_types();
  u = inverse_rosenblatt(u, num_threads);
  set_var_types_internal(actual_types);
  return u;
}

//! @brief calculates the log-likelihood.
//!
//! The log-likelihood is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, ..., U_{d, i}), \f]
//! where \f$ c \f$ is the copula density pdf().
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double
Vinecop::loglik(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  if (u.rows() < 1) {
    return this->get_loglik();
  } else {
    return pdf(u, num_threads).array().log().sum();
  }
}

//! @brief calculates the Akaike information criterion (AIC).
//!
//! The AIC is defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! get_npars(). The AIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double
Vinecop::aic(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  return -2 * this->loglik(u, num_threads) + 2 * get_npars();
}

//! @brief calculates the Bayesian information criterion (BIC).
//!
//! The BIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! get_npars(). The BIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double
Vinecop::bic(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  return -2 * this->loglik(u, num_threads) +
         get_npars() * log(static_cast<double>(u.rows()));
}

//! @brief calculates the modified Bayesian information criterion for vines
//! (mBICV).
//!
//! The mBICV is defined as
//! \f[ \mathrm{mBICV} = -2\, \mathrm{loglik} +  \ln(n) \nu, - 2 *
//! \sum_{t=1}^(d - 1) \{q_t log(\psi_0^t) - (d - t - q_t) log(1 -\psi_0^t)\}\f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood, \f$ \nu \f$ is the
//! (effective) number of parameters of the model, \f$ t \f$ is the tree level
//! \f$ \psi_0 \f$ is the prior probability of having a non-independence copula
//! in the first tree, and \f$ q_t \f$ is the number of non-independence copulas
//! in tree \f$ t \f$; The vBIC is a consistent model selection criterion for
//! parametric sparse vine copula models when \f$ d = o(\sqrt{n \ln n})\f$.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param psi0 baseline prior probability of a non-independence copula.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline double
Vinecop::mbicv(const Eigen::MatrixXd& u,
               const double psi0,
               const size_t num_threads) const
{

  double n = static_cast<double>(u.rows());
  double ll = this->loglik(u, num_threads);
  return -2 * ll + this->calculate_mbicv_penalty(n, psi0);
}

//! @brief returns sum of the number of parameters for all pair copulas (see
//! Bicop::get_npars()).
inline double
Vinecop::get_npars() const
{
  double npars = 0.0;
  for (auto& tree : pair_copulas_) {
    for (auto& pc : tree) {
      npars += pc.get_npars();
    }
  }
  return npars;
}

//! @brief calculates the Rosenblatt transform for a vine copula model.
//!
//! The Rosenblatt transform converts data from this model into independent
//! uniform variates. Only works for continuous data.
//!
//! @param data \f$ n \times d \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::MatrixXd
Vinecop::rosenblatt(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  if (get_n_discrete() > 0) {
    throw std::runtime_error("rosenblatt() only works for continuous models.");
  }
  check_data(u);
  size_t d = u.cols();
  size_t n = u.rows();

  // info about the vine structure (reverse rows (!) for more natural indexing)
  size_t trunc_lvl = vine_struct_.get_trunc_lvl();
  std::vector<size_t> order, inverse_order;
  TriangularArray<size_t> natural_array, min_array, needed_hfunc1,
    needed_hfunc2;
  if (trunc_lvl > 0) {
    order = vine_struct_.get_order();
    inverse_order = tools_stl::invert_permutation(order);
    natural_array = vine_struct_.get_struct_array(true);
    min_array = vine_struct_.get_min_array();
    needed_hfunc1 = vine_struct_.get_needed_hfunc1();
    needed_hfunc2 = vine_struct_.get_needed_hfunc2();
  }

  // fill first row of hfunc2 matrix with evaluation points;
  // points have to be reordered to correspond to natural order
  Eigen::MatrixXd hfunc1(n, d);
  Eigen::MatrixXd hfunc2(n, d);
  for (size_t j = 0; j < d; ++j)
    hfunc2.col(j) = u.col(order[j] - 1);

  auto do_batch = [&](const tools_batch::Batch& b) {
    Eigen::MatrixXd u_e(b.size, 2);
    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      tools_interface::check_user_interrupt(n * d > 1e5);
      for (size_t edge = 0; edge < d - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        size_t m = min_array(tree, edge);
        u_e.col(0) = hfunc2.block(b.begin, edge, b.size, 1);
        if (m == natural_array(tree, edge)) {
          u_e.col(1) = hfunc2.block(b.begin, m - 1, b.size, 1);
        } else {
          u_e.col(1) = hfunc1.block(b.begin, m - 1, b.size, 1);
        }

        // h-functions are only evaluated if needed in next step
        Bicop edge_copula = get_pair_copula(tree, edge).as_continuous();
        if (needed_hfunc1(tree, edge)) {
          hfunc1.block(b.begin, edge, b.size, 1) = edge_copula.hfunc1(u_e);
        }
        hfunc2.block(b.begin, edge, b.size, 1) = edge_copula.hfunc2(u_e);
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();
  }

  // go back to original order
  auto U_vine = u;
  for (size_t j = 0; j < d; j++) {
    U_vine.col(j) = hfunc2.col(inverse_order[j]);
  }

  return U_vine.array().min(1 - 1e-10).max(1e-10);
}

//! @brief calculates the inverse Rosenblatt transform for a vine copula model.
//!
//! The inverse Rosenblatt transform can be used for simulation: the
//! function applied to independent uniform variates resembles simulated
//! data from the vine copula model.
//!
//! If the problem is too large, it is split recursively into halves (w.r.t.
//! n, the number of observations).
//! "Too large" means that the required memory will exceed 1 GB. An
//! examplary configuration requiring less than 1 GB is \f$ n = 1000 \f$,
//! \f$d = 200\f$.
//!
//! Only works for continous models.
//!
//! @param u \f$ n \times d \f$ matrix of evaluation points.
//! @param num_threads the number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::MatrixXd
Vinecop::inverse_rosenblatt(const Eigen::MatrixXd& u,
                            const size_t num_threads) const
{
  if (get_n_discrete() > 0) {
    throw std::runtime_error(
      "inverse_rosenblatt() only works for continuous models.");
  }
  check_data(u);
  size_t n = u.rows();
  if (n < 1) {
    throw std::runtime_error("n must be at least one");
  }
  size_t d = d_;

  Eigen::MatrixXd U_vine = u.leftCols(d); // output matrix
  //                   (direct + indirect)    (U_vine)       (info matrices)
  size_t bytes_required = (8 * 2 * n * d * d) + (8 * n * d) + (4 * 4 * d * d);
  // if the problem is too large (requires more than 1 GB memory), split
  // the data into two halves and call simulate on the reduced data.
  if ((n > 1) & (bytes_required > 1e9)) {
    size_t n_half = n / 2;
    size_t n_left = n - n_half;
    U_vine.block(0, 0, n_half, d) =
      inverse_rosenblatt(u.block(0, 0, n_half, d));
    U_vine.block(n_half, 0, n_left, d) =
      inverse_rosenblatt(u.block(n_half, 0, n_left, d));
    return U_vine;
  }

  // info about the vine structure (in upper triangular matrix notation)
  size_t trunc_lvl = vine_struct_.get_trunc_lvl();
  std::vector<size_t> order, inverse_order;
  TriangularArray<size_t> natural_array, min_array, needed_hfunc1,
    needed_hfunc2;
  if (trunc_lvl > 0) {
    order = vine_struct_.get_order();
    inverse_order = tools_stl::invert_permutation(order);
    natural_array = vine_struct_.get_struct_array(true);
    min_array = vine_struct_.get_min_array();
    needed_hfunc1 = vine_struct_.get_needed_hfunc1();
    needed_hfunc2 = vine_struct_.get_needed_hfunc2();
  }

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects for (inverse) h-functions
    TriangularArray<Eigen::VectorXd> hinv2(d + 1, trunc_lvl + 1);
    TriangularArray<Eigen::VectorXd> hfunc1(d + 1, trunc_lvl + 1);

    // initialize with independent uniforms (corresponding to natural
    // order)
    for (size_t j = 0; j < d; ++j) {
      hinv2(std::min(trunc_lvl, d - j - 1), j) =
        u.block(b.begin, order[j] - 1, b.size, 1);
    }
    hfunc1(0, d - 1) = hinv2(0, d - 1);

    // loop through variables (0 is just the initial uniform)
    for (ptrdiff_t var = d - 2; var >= 0; --var) {
      tools_interface::check_user_interrupt(n * d > 1e5);
      size_t tree_start = std::min(trunc_lvl - 1, d - var - 2);
      for (ptrdiff_t tree = tree_start; tree >= 0; --tree) {
        Bicop edge_copula = get_pair_copula(tree, var).as_continuous();

        // extract data for conditional pair
        Eigen::MatrixXd U_e(b.size, 2);
        size_t m = min_array(tree, var);
        U_e.col(0) = hinv2(tree + 1, var);
        if (m == natural_array(tree, var)) {
          U_e.col(1) = hinv2(tree, m - 1);
        } else {
          U_e.col(1) = hfunc1(tree, m - 1);
        }

        // inverse Rosenblatt transform simulates data for conditional pair
        hinv2(tree, var) = edge_copula.hinv2(U_e);

        // if required at later stage, also calculate hfunc2
        if (var < static_cast<ptrdiff_t>(d_) - 1) {
          if (needed_hfunc1(tree, var)) {
            U_e.col(0) = hinv2(tree, var);
            hfunc1(tree + 1, var) = edge_copula.hfunc1(U_e);
          }
        }
      }
    }
    // go back to original order
    for (size_t j = 0; j < d; j++) {
      U_vine.block(b.begin, j, b.size, 1) = hinv2(0, inverse_order[j]);
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();
  }

  return U_vine;
}

//! checks if dimension d of the data matches the dimension of the vine.
inline void
Vinecop::check_data_dim(const Eigen::MatrixXd& data) const
{
  size_t d_data = data.cols();
  auto n_disc = get_n_discrete();
  size_t d_exp = d_ + n_disc;
  if ((d_data != d_exp) & (d_data != 2 * d_)) {
    std::stringstream msg;
    msg << "data has wrong number of columns; "
        << "expected: " << d_exp << " or " << 2 * d_ << ", actual: " << d_data
        << " (model contains ";
    if (n_disc == 0) {
      msg << "no discrete variables)." << std::endl;
    } else if (n_disc == 1) {
      msg << "1 discrete variable)." << std::endl;
    } else {
      msg << get_n_discrete() << "discrete variables)." << std::endl;
    }
    throw std::runtime_error(msg.str());
  }
}

//! checks if dimension d of the data matches the dimension of the vine.
inline void
Vinecop::check_data(const Eigen::MatrixXd& data) const
{
  check_data_dim(data);
  tools_eigen::check_if_in_unit_cube(data);
}

//! checks if pair copulas are compatible with the R-vine structure.
inline void
Vinecop::check_pair_copulas_rvine_structure(
  const std::vector<std::vector<Bicop>>& pair_copulas) const
{
  size_t trunc_lvl = vine_struct_.get_trunc_lvl();
  if (pair_copulas.size() > std::min(d_ - 1, trunc_lvl)) {
    std::stringstream message;
    message << "pair_copulas is too large; "
            << "expected size: < " << std::min(d_ - 1, trunc_lvl) << ", "
            << "actual size: " << pair_copulas.size() << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  for (size_t t = 0; t < pair_copulas.size(); ++t) {
    if (pair_copulas[t].size() != d_ - 1 - t) {
      std::stringstream message;
      message << "size of pair_copulas[" << t << "] "
              << "does not match dimension of matrix (" << d_ << "); "
              << "expected size: " << d_ - 1 - t << ", "
              << "actual size: " << pair_copulas[t].size() << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
Vinecop::finalize_fit(const tools_select::VinecopSelector& selector)
{
  vine_struct_ = selector.get_rvine_structure();
  threshold_ = selector.get_threshold();
  loglik_ = selector.get_loglik();
  nobs_ = selector.get_nobs();
  pair_copulas_ = selector.get_pair_copulas();
}

//! checks if weights are compatible with the data.
inline void
Vinecop::check_weights_size(const Eigen::VectorXd& weights,
                            const Eigen::MatrixXd& data) const
{
  if ((weights.size() > 0) & (weights.size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! checks if data size is large enough
inline void
Vinecop::check_enough_data(const Eigen::MatrixXd& data) const
{
  if (data.rows() == 1) {
    throw std::runtime_error("data must have more than one row");
  }
}

inline void
Vinecop::check_fitted() const
{
  if (std::isnan(loglik_)) {
    throw std::runtime_error("copula has not been fitted from data ");
  }
}

//! @brief truncate the vine copula model.
//! @param trunc_lvl the truncation level.
//! If the model is already truncated at a level less than `trunc_lvl`,
//! the function does nothing.
inline void
Vinecop::truncate(size_t trunc_lvl)
{
  if (trunc_lvl < this->get_trunc_lvl()) {
    vine_struct_.truncate(trunc_lvl);
    pair_copulas_.resize(trunc_lvl);
  }
}

//! set all variable types to continuous.
//! the function can be const, because var_types_ is mutable.
inline void
Vinecop::set_continuous_var_types() const
{
  var_types_ = std::vector<std::string>(d_);
  for (auto& t : var_types_)
    t = "c";
  set_var_types_internal(var_types_);
}

//! returns the number of discrete variables.
inline int
Vinecop::get_n_discrete() const
{
  int n_discrete = 0;
  for (auto t : var_types_) {
    n_discrete += (t == "d");
  }
  return n_discrete;
}

//! removes superfluous columns for continuous data.
inline Eigen::MatrixXd
Vinecop::collapse_data(const Eigen::MatrixXd& u) const
{
  if (static_cast<size_t>(u.cols()) == d_ + get_n_discrete()) {
    return u;
  }
  Eigen::MatrixXd u_new(u.rows(), d_ + get_n_discrete());
  u_new.leftCols(d_) = u.leftCols(d_);
  size_t disc_count = 0;
  for (size_t i = 0; i < d_; ++i) {
    if (var_types_[i] == "d") {
      u_new.col(d_ + disc_count++) = u.col(d_ + i);
    }
  }
  return u_new;
}

//! summarizes the model into a string (can be used for printing).
inline std::string
Vinecop::str() const
{
  std::stringstream str;
  auto arr = vine_struct_.get_struct_array();
  auto order = vine_struct_.get_order();
  for (size_t t = 0; t < vine_struct_.get_trunc_lvl(); ++t) {
    str << "** Tree: " << t << std::endl;
    for (size_t e = 0; e < d_ - 1 - t; ++e) {
      str << order[e] << "," << arr(t, e);
      if (t > 0) {
        str << " | ";
        for (size_t cv = t - 1; cv > 1; --cv) {
          str << arr(cv, e) - 1 << ",";
        }
        str << arr(0, e);
      }
      str << " <-> " << pair_copulas_[t][e].str() << std::endl;
    }
  }
  return str.str();
}
}
