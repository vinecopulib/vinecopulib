// Copyright © 2016-2020 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/bicop/indep.hpp>

namespace vinecopulib {
//! @brief An abstract class for rotation mixtures of copula families.
//!
//! For a given base model \f$ C(\cdot; \theta) \$f with parameter $\theta$, the
//! rotation mixture family is given by
//! \f$ (C_0(\cdot; \theta_0) + C_{90}(\cdot; \theta_{90}) + C_{180}(\cdot;
//! \theta_{180}) + C_{270}(\cdot; \theta_{270})) / 4 \f$, where \f$ C_0, C_{90},
//! C_{180},\f$ and C_{270} are rotated versions of $C$. 
//! 
//! This class is used in the implementation underlying the Bicop class. Users 
//! should not use RotationMixtureBicop or derived classes directly, but always 
//! work with the Bicop interface.
class RotationMixtureBicop : public ParBicop
{
public:
  // constructor
  RotationMixtureBicop();

protected:
  BicopPtr base_bicop_;

private:
  // PDF
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u);

  // CDF
  Eigen::VectorXd cdf(const Eigen::MatrixXd& u);

  // hfunction
  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u);
  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u);

  // inverse hfunction
  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u);
  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u);

  Eigen::MatrixXd tau_to_parameters(const double& tau);
  double parameters_to_tau(const Eigen::MatrixXd& parameters);

  void rotate_90(Eigen::MatrixXd& u) const;

};

//! @brief Rotation mixture of Gumbel copulas.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class RMGumbelBicop : public RotationMixtureBicop
{
public:
  // constructor
  RMGumbelBicop();

private:
  Eigen::VectorXd get_start_parameters(const double tau);
};


//! @brief Rotation mixture of Clayton copulas.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class RMClaytonBicop : public RotationMixtureBicop
{
public:
  // constructor
  RMClaytonBicop();

private:
  Eigen::VectorXd get_start_parameters(const double tau);
};


//! @brief Rotation mixture of Joe copulas.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class RMJoeBicop : public RotationMixtureBicop
{
public:
  // constructor
  RMJoeBicop();

private:
  Eigen::VectorXd get_start_parameters(const double tau);
};

}

#include <vinecopulib/bicop/implementation/rotation_mixture.ipp>
