// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/fit_controls.hpp>
#include <vinecopulib/misc/nlohmann_json.hpp>

namespace vinecopulib {

// forward declaration of Abstract class
class AbstractBicop;
using BicopPtr = std::shared_ptr<AbstractBicop>;

//! @brief A class for bivariate copula models.
//!
//! The copula model is fully characterized by the family, rotation,
//! and parameters.
class Bicop
{

public:
  // Constructors
  Bicop(const BicopFamily family = BicopFamily::indep,
        const int rotation = 0,
        const Matrix& parameters = Matrix(),
        const std::vector<std::string>& var_types = { "c", "c" });

  explicit Bicop(const Matrix& data,
                 const FitControlsBicop& controls = FitControlsBicop(),
                 const std::vector<std::string>& var_types = { "c", "c" });

  Bicop(const Bicop& other);

  explicit Bicop(const std::string& filename);

  explicit Bicop(const nlohmann::json& input);

  Bicop& operator=(Bicop other);

  // Serialize
  nlohmann::json to_json() const;

  void to_file(const std::string& filename) const;

  // Getters and setters
  BicopFamily get_family() const;

  std::string get_family_name() const;

  int get_rotation() const;

  Matrix get_parameters() const;

  double get_tau() const;

  double get_npars() const;

  double get_loglik() const;
  size_t get_nobs() const;
  double get_aic() const;
  double get_bic() const;
  double get_mbic(const double psi0 = 0.9) const;

  void set_rotation(const int rotation);

  void set_parameters(const Matrix& parameters);

  void set_var_types(const std::vector<std::string>& var_types = { "c", "c" });

  std::vector<std::string> get_var_types() const;

  // Stats methods
  Vector pdf(const Matrix& u) const;

  Vector cdf(const Matrix& u) const;

  Vector hfunc1(const Matrix& u) const;

  Vector hfunc2(const Matrix& u) const;

  Vector hinv1(const Matrix& u) const;

  Vector hinv2(const Matrix& u) const;

  Matrix simulate(const size_t& n,
                  const bool qrng = false,
                  const std::vector<int>& seeds = std::vector<int>()) const;

  // Methods modifying the family/rotation/parameters
  void fit(const Matrix& data,
           const FitControlsBicop& controls = FitControlsBicop());

  void select(const Matrix& data,
              FitControlsBicop controls = FitControlsBicop());

  // Fit statistics
  double loglik(const Matrix& u = Matrix()) const;

  double aic(const Matrix& u = Matrix()) const;

  double bic(const Matrix& u = Matrix()) const;

  double mbic(const Matrix& u = Matrix(), const double psi0 = 0.9) const;

  // Misc
  std::string str() const;

  double parameters_to_tau(const Matrix& parameters) const;

  Matrix tau_to_parameters(const double& tau) const;

  void flip();

  Matrix get_parameters_lower_bounds() const;

  Matrix get_parameters_upper_bounds() const;

  Bicop as_continuous() const;

private:
  Matrix format_data(const Matrix& u) const;

  void rotate_data(Matrix& u) const;

  Matrix prep_for_abstract(const Matrix& u) const;

  void check_rotation(int rotation) const;

  void check_data(const Matrix& u) const;

  void check_data_dim(const Matrix& u) const;

  void check_var_types(const std::vector<std::string>& var_types) const;

  void flip_abstract_var_types();

  void check_weights_size(const Vector& weights,
                          const Matrix& data) const;

  void check_fitted() const;

  unsigned short get_n_discrete() const;

  double compute_mbic_penalty(const size_t nobs, const double psi0) const;

  BicopPtr get_bicop() const;

  BicopPtr bicop_;
  int rotation_{ 0 };
  size_t nobs_{ 0 };
  mutable std::vector<std::string> var_types_;
};
}

#include <vinecopulib/bicop/implementation/class.ipp>
