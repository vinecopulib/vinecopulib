// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <memory>

#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_linalg.hpp>

namespace vinecopulib {
//! @brief An abstract class for bivariate copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class AbstractBicop
{
  friend class Bicop;

public:
  virtual ~AbstractBicop() = 0;

protected:
  // Factories
  static std::shared_ptr<AbstractBicop> create(
    BicopFamily family = BicopFamily::indep,
    const Matrix& parameters = Matrix());

  // Getters and setters
  BicopFamily get_family() const;

  std::string get_family_name() const;

  double get_loglik() const;

  void set_loglik(const double loglik = NAN);

  void set_var_types(const std::vector<std::string>& var_types);

  virtual Matrix get_parameters() const = 0;

  virtual Matrix get_parameters_lower_bounds() const = 0;

  virtual Matrix get_parameters_upper_bounds() const = 0;

  virtual void set_parameters(const Matrix& parameters) = 0;

  // Virtual methods
  virtual void fit(const Matrix& data,
                   std::string method,
                   double mult,
                   const Eigen::VectorXd& weights) = 0;

  virtual double get_npars() const = 0;

  virtual void set_npars(const double& npars) = 0;

  virtual double parameters_to_tau(const Matrix& parameters) = 0;

  virtual void flip() = 0;

  // following are virtual so they can be overriden by KernelBicop
  virtual Eigen::VectorXd pdf(const Matrix& u);

  virtual Eigen::VectorXd cdf(const Matrix& u) = 0;

  virtual Eigen::VectorXd hfunc1(const Matrix& u);

  virtual Eigen::VectorXd hfunc2(const Matrix& u);

  Eigen::VectorXd hinv1(const Matrix& u);

  Eigen::VectorXd hinv2(const Matrix& u);

  virtual Eigen::VectorXd pdf_raw(const Matrix& u) = 0;

  virtual Eigen::VectorXd hfunc1_raw(const Matrix& u) = 0;

  virtual Eigen::VectorXd hfunc2_raw(const Matrix& u) = 0;

  virtual Eigen::VectorXd hinv1_raw(const Matrix& u) = 0;

  virtual Eigen::VectorXd hinv2_raw(const Matrix& u) = 0;

  virtual Matrix tau_to_parameters(const double& tau) = 0;
  Matrix no_tau_to_parameters(const double&);

  // Misc methods
  Eigen::VectorXd hinv1_num(const Matrix& u);

  Eigen::VectorXd hinv2_num(const Matrix& u);

  Eigen::VectorXd pdf_c_d(const Matrix& u);

  Eigen::VectorXd pdf_d_d(const Matrix& u);

  double loglik(const Matrix& u,
                const Eigen::VectorXd weights = Eigen::VectorXd());

  // Data members
  BicopFamily family_;
  double loglik_{ NAN };
  std::vector<std::string> var_types_{ "c", "c" };
};

//! A shared pointer to an object of class AbstracBicop.
typedef std::shared_ptr<AbstractBicop> BicopPtr;
}

#include <vinecopulib/bicop/implementation/abstract.ipp>
