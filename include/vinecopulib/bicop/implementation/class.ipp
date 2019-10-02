// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <mutex>
#include <vinecopulib/bicop/abstract.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {

//! @brief creates a specific bivariate copula model.
//! @param family the copula family.
//! @param rotation the rotation of the copula; one of 0, 90, 180, or 270
//!     (for Independence, Gaussian, Student, Frank, and nonparametric
//!     families, only 0 is allowed).
//! @param parameters the copula parameters.
inline Bicop::Bicop(const BicopFamily family,
                    const int rotation,
                    const Eigen::MatrixXd& parameters)
{
  bicop_ = AbstractBicop::create(family, parameters);
  // family must be set before checking the rotation
  set_rotation(rotation);
  if (bicop_->get_family() != BicopFamily::indep) {
    bicop_->set_loglik();
  } else {
    bicop_->set_loglik(0.0);
  }
}

//! @brief create a copula model from the data,
//! equivalent to `Bicop cop; cop.select(data, controls)`.
//!
//! @param data see select().
//! @param controls see select().
inline Bicop::Bicop(const Eigen::MatrixXd& data,
                    const FitControlsBicop& controls)
{
  select(data, controls);
}

//! @brief creates from a boost::property_tree::ptree object
//! @param input the boost::property_tree::ptree object to convert from
//! (see to_ptree() for the structure of the input).
inline Bicop::Bicop(const boost::property_tree::ptree input)
  : Bicop(get_family_enum(input.get<std::string>("family")),
          input.get<int>("rotation"),
          tools_serialization::ptree_to_matrix<double>(
            input.get_child("parameters")))
{
  var_types_ = tools_serialization::ptree_to_vector<std::string>(
    input.get_child("var_types"));
  // try block for backwards compatibility
  try {
    nobs_ = input.get<size_t>("nobs_");
    bicop_->set_loglik(input.get<double>("loglik"));
  } catch (...) {}
}

//! @brief creates from a JSON file
//! @param filename the name of the JSON file to read (see to_ptree() for the
//! structure of the file).
inline Bicop::Bicop(const std::string filename)
  : Bicop(tools_serialization::json_to_ptree(filename.c_str()))
{}

//! @brief Convert the copula into a boost::property_tree::ptree object
//!
//! The boost::property_tree::ptree is contains of three values named
//! `"family"`, `"rotation"`, `"parameters"`, respectively a string
//! for the family name, an integer for the rotation, and an Eigen::MatrixXd
//! for the parameters.
//!
//! @return the boost::property_tree::ptree object containing the copula.
inline boost::property_tree::ptree
Bicop::to_ptree() const
{
  boost::property_tree::ptree output;

  output.put("family", get_family_name());
  output.put("rotation", rotation_);
  auto mat_node = tools_serialization::matrix_to_ptree(get_parameters());
  output.add_child("parameters", mat_node);
  output.add_child("var_types",
                   tools_serialization::vector_to_ptree(var_types_));

  output.put("nobs_", nobs_);
  output.put("loglik", bicop_->get_loglik());

  return output;
}

//! @brief Write the copula object into a JSON file
//!
//! See to_ptree() for the structure of the file.
//!
//! @param filename the name of the file to write.
inline void
Bicop::to_json(const std::string filename) const
{
  boost::property_tree::write_json(filename.c_str(), to_ptree());
}

//! @brief evaluates the copula density.
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The copula density evaluated at \c u.
inline Eigen::VectorXd
Bicop::pdf(const Eigen::MatrixXd& u) const
{
  check_data(u);
  return bicop_->pdf(prep_for_abstract(u));
}

//! @brief evaluates the copula distribution.
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The copula distribution evaluated at \c u.
inline Eigen::VectorXd
Bicop::cdf(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd p = bicop_->cdf(prep_for_abstract(u).leftCols(2));
  switch (rotation_) {
    default:
      return p;

    case 90:
      return u.col(1) - p;

    case 180:
      return p.array() - 1 + u.leftCols(2).rowwise().sum().array();

    case 270:
      return u.col(0) - p;
  }
}

//! @brief calculates the first h-function.
//!
//! The first h-function is
//! \f$ h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1) \f$.
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hfunc1(const Eigen::MatrixXd& u) const
{
  check_data(u);
  switch (rotation_) {
    default:
      return bicop_->hfunc1(prep_for_abstract(u));

    case 90:
      return bicop_->hfunc2(prep_for_abstract(u));

    case 180:
      return 1.0 - bicop_->hfunc1(prep_for_abstract(u)).array();

    case 270:
      return 1.0 - bicop_->hfunc2(prep_for_abstract(u)).array();
  }
}

//! @brief calculates the second h-function.
//!
//! The second h-function is
//! \f$ h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)  \f$.
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hfunc2(const Eigen::MatrixXd& u) const
{
  check_data(u);
  switch (rotation_) {
    default:
      return bicop_->hfunc2(prep_for_abstract(u));

    case 90:
      return 1.0 - bicop_->hfunc1(prep_for_abstract(u)).array();

    case 180:
      return 1.0 - bicop_->hfunc2(prep_for_abstract(u)).array();

    case 270:
      return bicop_->hfunc1(prep_for_abstract(u));
  }
}

//! @brief calculates the inverse of \f$ h_1 \f$ (see hfunc1()) w.r.t. the
//! second argument.
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hinv1(const Eigen::MatrixXd& u) const
{
  check_data(u);
  switch (rotation_) {
    default:
      return bicop_->hinv1(prep_for_abstract(u));

    case 90:
      return bicop_->hinv2(prep_for_abstract(u));

    case 180:
      return 1.0 - bicop_->hinv1(prep_for_abstract(u)).array();

    case 270:
      return 1.0 - bicop_->hinv2(prep_for_abstract(u)).array();
  }
}

//! @brief calculates the inverse of \f$ h_2 \f$ (see hfunc2()) w.r.t. the first
//! argument.
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hinv2(const Eigen::MatrixXd& u) const
{
  check_data(u);
  switch (rotation_) {
    default:
      return bicop_->hinv2(prep_for_abstract(u));

    case 90:
      return 1.0 - bicop_->hinv1(prep_for_abstract(u)).array();

    case 180:
      return 1.0 - bicop_->hinv2(prep_for_abstract(u)).array();

    case 270:
      return bicop_->hinv1(prep_for_abstract(u));
  }
}
//! @}

//! @brief simulates from a bivariate copula.
//!
//! @param n number of observations.
//! @param qrng set to true for quasi-random numbers.
//! @param seeds seeds of the (quasi-)random number generator; if empty
//! (default),
//!   the (quasi-)random number generator is seeded randomly.
//! @return An \f$ n \times 2 \f$ matrix of samples from the copula model.
inline Eigen::MatrixXd
Bicop::simulate(const size_t& n,
                const bool qrng,
                const std::vector<int>& seeds) const
{
  auto u = tools_stats::simulate_uniform(n, 2, qrng, seeds);
  // use inverse Rosenblatt transform to generate a sample from the copula
  // (always simulate continuous data)
  u.col(1) = this->as_continuous().hinv1(u);
  return u;
}

//! @brief calculates the log-likelihood.
//!
//! The log-likelihood is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, U_{2, i}), \f]
//! where \f$ c \f$ is the copula density pdf().
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline double
Bicop::loglik(const Eigen::MatrixXd& u) const
{
  if (u.rows() < 1) {
    return get_loglik();
  } else {
    tools_eigen::check_if_in_unit_cube(u);
    return bicop_->loglik(prep_for_abstract(u));
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
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline double
Bicop::aic(const Eigen::MatrixXd& u) const
{
  return -2 * loglik(u) + 2 * get_npars();
}

//! @brief calculates the Bayesian information criterion (BIC).
//!
//! The BIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! get_npars(). The BIC is a consistent model selection criterion
//! for parametric models.
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
inline double
Bicop::bic(const Eigen::MatrixXd& u) const
{
  Eigen::MatrixXd u_no_nan = u;
  double n = nobs_;
  if (u.rows() > 0) {
    tools_eigen::remove_nans(u_no_nan);
    n = static_cast<double>(u_no_nan.rows());
  }
  return -2 * loglik(u_no_nan) + get_npars() * log(n);
}

//! @brief calculates the modified Bayesian information criterion (mBIC).
//!
//! The mBIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \nu \ln(n)
//!  - 2 (I log(\psi_0) + (1 - I) log(1 - \psi_0) \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ \nu \f$ is the
//! (effective) number of parameters of the model, \f$ \psi_0 \f$ is the prior
//! probability of having a non-independence copula and \f$ I \f$ is an
//! indicator for the family being non-independence; see loglik() and
//! get_npars().
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
//! @param psi0 prior probability of a non-independence copula.
inline double
Bicop::mbic(const Eigen::MatrixXd& u, const double psi0) const
{
  Eigen::MatrixXd u_no_nan = u;
  bool is_indep = (this->get_family() == BicopFamily::indep);
  double npars = this->get_npars();
  double log_prior = static_cast<double>(!is_indep) * std::log(psi0) +
                     static_cast<double>(is_indep) * std::log(1.0 - psi0);
  double n = nobs_;
  if (u.rows() > 0) {
    n = static_cast<double>(u_no_nan.rows());
  }
  return -2 * this->loglik(u_no_nan) + std::log(n) * npars - 2 * log_prior;
}

//! @brief returns the actual number of parameters for parameteric families.
//!
//! For nonparametric families, there is a conceptually similar definition in
//! the sense that it can be used in the calculation of fit statistics.
inline double
Bicop::get_npars() const
{
  return bicop_->get_npars();
}

//! @brief converts a Kendall's \f$ \tau \f$ to the copula parameters of the
//! current family
//!
//! (only works for one-parameter families)
//! @param tau a value in \f$ (-1, 1) \f$.
inline Eigen::MatrixXd
Bicop::tau_to_parameters(const double& tau) const
{
  return bicop_->tau_to_parameters(tau);
}

//! @brief converts the parameters to the Kendall's \f$ tau \f$ for the current
//! family.
//!
//! @param parameters the parameters (must be a valid parametrization of
//!     the current family).
inline double
Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters) const
{
  double tau = bicop_->parameters_to_tau(parameters);
  if (tools_stl::is_member(rotation_, { 90, 270 })) {
    tau *= -1;
  }
  return tau;
}

//! @name Getters and setters
//!
//! @{

//! get the copula family
inline BicopFamily
Bicop::get_family() const
{
  return bicop_->get_family();
}

//! get the copula family as a string
inline std::string
Bicop::get_family_name() const
{
  return bicop_->get_family_name();
}

//! get the rotation
inline int
Bicop::get_rotation() const
{
  return rotation_;
}

//! get the parameters
inline Eigen::MatrixXd
Bicop::get_parameters() const
{
  return bicop_->get_parameters();
}

//! get the log-likelihood (only for fitted objects)
inline double
Bicop::get_loglik() const
{
  check_fitted();
  return bicop_->get_loglik();
}

//! get the number of observations (only for fitted objects)
inline size_t
Bicop::get_nobs() const
{
  check_fitted();
  return nobs_;
}

//! get the aic (only for fitted objects)
inline double
Bicop::get_aic() const
{
  check_fitted();
  return -2 * bicop_->get_loglik() + 2 * bicop_->get_npars();
}

//! get the bic (only for fitted objects)
inline double
Bicop::get_bic() const
{
  check_fitted();
  double npars = bicop_->get_npars();
  return -2 * bicop_->get_loglik() + std::log(nobs_) * npars;
}

//! get the modified bic (only for fitted objects)
inline double
Bicop::get_mbic(const double psi0) const
{
  check_fitted();
  return -2 * bicop_->get_loglik() + compute_mbic_penalty(nobs_, psi0);
}

inline double
Bicop::compute_mbic_penalty(const size_t nobs, const double psi0) const
{
  double npars = bicop_->get_npars();
  bool is_indep = (this->get_family() == BicopFamily::indep);
  double log_prior = static_cast<double>(!is_indep) * std::log(psi0) +
                     static_cast<double>(is_indep) * std::log(1.0 - psi0);
  return std::log(nobs) * npars - 2 * log_prior;
}

//! get the Kendall's tau
inline double
Bicop::get_tau() const
{
  return parameters_to_tau(bicop_->get_parameters());
}

//! set the rotation
inline void
Bicop::set_rotation(const int rotation)
{
  check_rotation(rotation);
  if ((rotation_ - rotation % 180) != 0) {
    flip_var_types();
  }
  rotation_ = rotation;
  bicop_->set_loglik();
}

inline void
Bicop::check_data(const Eigen::MatrixXd& u) const
{
  check_data_dim(u);
  tools_eigen::check_if_in_unit_cube(u);
}

inline void
Bicop::check_data_dim(const Eigen::MatrixXd& u) const
{
  size_t n_cols = u.cols();
  size_t n_cols_exp = 2 + get_n_discrete();
  if ((n_cols != n_cols_exp) & (n_cols != 4)) {
    std::stringstream msg;
    msg << "data has wrong number of columns; "
        << "expected: " << n_cols_exp << " or 4, actual: " << n_cols 
        << " (model contains ";
    if (n_cols_exp == 2) {
      msg << "no ";
    } else {
      msg << get_n_discrete() << " ";
    }
    msg << "discrete variables)." << std::endl;
    throw std::runtime_error(msg.str());
  }
}

inline void
Bicop::flip_var_types()
{
  std::swap(bicop_->var_types_[0], bicop_->var_types_[1]);
}

inline void
Bicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  bicop_->set_parameters(parameters);
  bicop_->set_loglik();
}

//! @brief sets variable types.
//! @param var_types a vector of size two specifying the types of the variables,
//!   e.g., `{"c", "d"}` means first varible continuous, second discrete.
inline void
Bicop::set_var_types(const std::vector<std::string>& var_types)
{
  if (var_types.size() != 2) {
    throw std::runtime_error("var_types must have size two.");
  }
  for (auto t : var_types) {
    if (!tools_stl::is_member(t, { "c", "d" })) {
      throw std::runtime_error("var type must be either 'c' or 'd'.");
    }
  }
  var_types_ = var_types;
  bicop_->set_var_types(var_types);
  if (tools_stl::is_member(static_cast<size_t>(rotation_), { 90, 270 })) {
    flip_var_types();
  }
}

//! @brief extracts variable types.
inline std::vector<std::string>
Bicop::get_var_types() const
{
  return var_types_;
}
//! @}

//! @name Utilities
//! @{
//! useful functions for bivariate copulas

//! adjust's the copula model to a change in the variable order.
inline void
Bicop::flip()
{
  BicopFamily family = bicop_->get_family();
  if (tools_stl::is_member(family, bicop_families::flip_by_rotation)) {
    double loglik = bicop_->get_loglik();
    if (rotation_ == 90) {
      set_rotation(270);
    } else if (rotation_ == 270) {
      set_rotation(90);
    }
    bicop_->set_loglik(loglik);
  } else {
    flip_var_types();
    bicop_->flip();
  }
}

//! summarizes the model into a string (can be used for printing).
inline std::string
Bicop::str() const
{
  std::stringstream bicop_str;
  bicop_str << get_family_name();
  if (get_rotation() != 0) {
    bicop_str << " " << get_rotation() << "°";
  }
  if (get_family() == BicopFamily::tll) {
    bicop_str << ", parameters = [30x30 grid]";
  } else if (get_family() != BicopFamily::indep) {
    bicop_str << ", parameters = " << get_parameters();
  }
  return bicop_str.str().c_str();
}

//! extract lower bounds for copula parameters.
inline Eigen::MatrixXd
Bicop::get_parameters_lower_bounds() const
{
  return bicop_->get_parameters_lower_bounds();
}

//! extract upper bounds for copula parameters.
inline Eigen::MatrixXd
Bicop::get_parameters_upper_bounds() const
{
  return bicop_->get_parameters_upper_bounds();
}

//! @}

inline BicopPtr
Bicop::get_bicop() const
{
  return bicop_;
}

inline Bicop
Bicop::as_continuous() const
{
  std::vector<std::string> cc = { "c", "c" };
  if (var_types_ == cc)
    return *this;
  auto bc_new = *this;
  bc_new.set_var_types(cc);
  return bc_new;
}

//! fits a bivariate copula (with fixed family) to data.
//!
//! For parametric models, two different methods are available. `"mle"` fits
//! the parameters by maximum-likelihood. `"itau"` uses inversion of
//! Kendall's \f$ \tau \f$, but is only available for one-parameter families
//! and the Student t copula. For the latter, there is a one-to-one
//! transformation for the first parameter, the second is found by profile
//! likelihood optimization (with accuracy of at least 0.5). Nonparametric
//! families have specialized methods, no specification is required.
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ F_{X_1}(X_1), F_{X_2}(X_2) \f$. Let \f$ k \f$ denote the number of
//! discrete variables (either one or two). Then the second \f$ n \times k \f$
//! block contains realizations of \f$ F_{X_k}(X_k^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_{X_k}(X_k^-) = F_{X_k}(X_k - 1) \f$.
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
//! @param controls the controls (see FitControlsBicop).
inline void
Bicop::fit(const Eigen::MatrixXd& data, const FitControlsBicop& controls)
{
  std::string method;
  if (tools_stl::is_member(bicop_->get_family(), bicop_families::parametric)) {
    method = controls.get_parametric_method();
  } else {
    method = controls.get_nonparametric_method();
  }
  tools_eigen::check_if_in_unit_cube(data);

  auto w = controls.get_weights();
  Eigen::MatrixXd data_no_nan = data;
  check_weights_size(w, data);
  tools_eigen::remove_nans(data_no_nan, w);

  bicop_->fit(prep_for_abstract(data_no_nan),
              method,
              controls.get_nonparametric_mult(),
              w);
  nobs_ = data_no_nan.rows();
}

//

//! @brief selects the best fitting model.
//!
//! The function calls fit() for all families in
//! `family_set` and selecting the best fitting model by either BIC or AIC,
//! see bic() and aic().
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ F_{X_1}(X_1), F_{X_2}(X_2) \f$. Let \f$ k \f$ denote the number of
//! discrete variables (either one or two). Then the second \f$ n \times k \f$
//! block contains realizations of \f$ F_{X_k}(X_k^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_{X_k}(X_k^-) = F_{X_k}(X_k - 1) \f$.
//!
//! @param data an \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1)^2 \f$, where \f$ k \f$ is the number of discrete variables.
//! @param controls the controls (see FitControlsBicop).
inline void
Bicop::select(const Eigen::MatrixXd& data, FitControlsBicop controls)
{
  using namespace tools_select;
  check_weights_size(controls.get_weights(), data);
  Eigen::MatrixXd data_no_nan = data;
  {
    auto w = controls.get_weights();
    tools_eigen::remove_nans(data_no_nan, w);
    controls.set_weights(w);
  }
  check_data(data_no_nan);
  nobs_ = data_no_nan.rows();

  bicop_ = AbstractBicop::create();
  bicop_->set_var_types(var_types_);
  rotation_ = 0;
  bicop_->set_loglik(0.0);
  if (data_no_nan.rows() >= 10) {
    data_no_nan = clip_data(data_no_nan);
    std::vector<Bicop> bicops = create_candidate_bicops(data_no_nan, controls);
    for (auto& bc : bicops) {
      bc.set_var_types(var_types_);
    }

    // Estimate all models and select the best one using the
    // selection_criterion
    double fitted_criterion = std::numeric_limits<double>::max();
    std::mutex m;
    auto fit_and_compare = [&](Bicop cop) {
      tools_interface::check_user_interrupt();
      // Estimate the model
      cop.fit(data_no_nan, controls);

      // Compute the selection criterion
      double new_criterion;
      double ll = cop.get_loglik();
      if (controls.get_selection_criterion() == "loglik") {
        new_criterion = -ll;
      } else if (controls.get_selection_criterion() == "aic") {
        new_criterion = -2 * ll + 2 * cop.get_npars();
      } else {
        double n_eff = static_cast<double>(data_no_nan.rows());
        if (controls.get_weights().size() > 0) {
          n_eff = std::pow(controls.get_weights().sum(), 2);
          n_eff /= controls.get_weights().array().pow(2).sum();
        }
        double npars = cop.get_npars();

        new_criterion = -2 * ll + log(n_eff) * npars; // BIC
        if (controls.get_selection_criterion() == "mbic") {
          // correction for mBIC
          bool is_indep = (this->get_family() == BicopFamily::indep);
          double psi0 = controls.get_psi0();
          double log_prior = static_cast<double>(!is_indep) * log(psi0) +
                             static_cast<double>(is_indep) * log(1.0 - psi0);
          new_criterion -= 2 * log_prior;
        }
      }

      // the following block modifies thread-external variables
      // and is thus shielded by a mutex
      {
        std::lock_guard<std::mutex> lk(m);
        // If the new model is better than the current one,
        // then replace the current model by the new one
        if (new_criterion < fitted_criterion) {
          fitted_criterion = new_criterion;
          bicop_ = cop.get_bicop();
          rotation_ = cop.get_rotation();
        }
      }
    };

    tools_thread::ThreadPool pool(controls.get_num_threads());
    pool.map(fit_and_compare, bicops);
  }
}

//! adds an additional column if there's only one discrete variable. 
//! (continuous models only require two columns, discrete models always four)
inline Eigen::MatrixXd
Bicop::extend_data(const Eigen::MatrixXd& u) const
{
  if ((get_n_discrete() != 1) | (u.cols() == 4)) {
    return u;
  }
  Eigen::MatrixXd u_new(u.rows(), 4);
  u_new.leftCols(2) = u.leftCols(2);
  int disc_col = (var_types_[1] == "d");
  int cont_col = 1 - disc_col;
  u_new.col(2 + disc_col) = u.col(2);
  u_new.col(2 + cont_col) = u.col(cont_col);
  return u_new;
}

//! clisp the data to the interval [1e-10, 1 - 1e-10] for numerical stability.
inline Eigen::MatrixXd
Bicop::clip_data(const Eigen::MatrixXd& u) const
{
  auto clip = [] (const double& x) {
    return std::min(std::max(x, 1e-10), 1 - 1e-10);
  };
  return tools_eigen::unaryExpr_or_nan(u, clip);
}

//! rotates the data corresponding to the models rotation.
//! @param u an `n x 2` matrix.
inline Eigen::MatrixXd
Bicop::rotate_data(const Eigen::MatrixXd& u) const
{
  auto u_new = u;
  // counter-clockwise rotations
  switch (rotation_) {
    case 0:
      break;

    case 90:
      u_new.col(0) = u.col(1);
      u_new.col(1) = 1.0 - u.col(0).array();
      break;

    case 180:
      u_new.col(0) = 1.0 - u.col(0).array();
      u_new.col(1) = 1.0 - u.col(1).array();
      break;

    case 270:
      u_new.col(0) = 1.0 - u.col(1).array();
      u_new.col(1) = u.col(0);
      break;
  }
  return u_new;
}

//! prepares data for use with the `AbstractBicop` class:
//! - add an additional column if there's only one discrete variable.
//! - clip the data to the interval [1e-10, 1 - 1e-10] for numerical stability.
//! - rotate the data appropriately (`AbstractBicop` is always 0deg-rotation).
inline Eigen::MatrixXd
Bicop::prep_for_abstract(const Eigen::MatrixXd& u) const
{
  auto u_new = extend_data(u);
  u_new = clip_data(u_new);
  u_new.leftCols(2) = rotate_data(u_new.leftCols(2));
  if (u_new.cols() > 2) {
    u_new.rightCols(2) = rotate_data(u_new.rightCols(2));
  }
  return u_new;
}

//! checks whether the supplied rotation is valid (only 0, 90, 180, 270 allowd).
inline void
Bicop::check_rotation(int rotation) const
{
  using namespace tools_stl;
  std::vector<int> allowed_rotations = { 0, 90, 180, 270 };
  if (!is_member(rotation, allowed_rotations)) {
    throw std::runtime_error("rotation must be one of {0, 90, 180, 270}");
  }
  if (is_member(bicop_->get_family(), bicop_families::rotationless)) {
    if (rotation != 0) {
      throw std::runtime_error("rotation must be 0 for the " +
                               bicop_->get_family_name() + " copula");
    }
  }
}

//! checks whether weights and data have matching sizes.
inline void
Bicop::check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const
{
  if ((weights.size() > 0) & (weights.size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! checks whether the Bicop object was fitted to data.
inline void
Bicop::check_fitted() const
{
  if ((boost::math::isnan)(bicop_->get_loglik())) {
    throw std::runtime_error("copula has not been fitted from data or its "
                             "parameters have been modified manually");
  }
}

//! returns the number of discrete variables.
inline int
Bicop::get_n_discrete() const
{
  int n_discrete = 0;
  for (auto t : var_types_) {
    n_discrete += (t == "d");
  }
  return n_discrete;
}

}
