// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

inline VonMisesBicop::VonMisesBicop()
{
  family_ = BicopFamily::von_mises;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  const double inf = std::numeric_limits<double>::infinity();
  parameters_ << 0, 0;
  parameters_lower_bounds_ << 0, -inf;
  parameters_upper_bounds_ << 100, inf;
}

inline const std::vector<double>&
VonMisesBicop::ratios_for(Series& series, double kappa)
{
  if (series.kappa != kappa) {
    series.kappa = kappa;
    series.ratios = tools_circular::bessel_i_ratios(
      kappa, tools_circular::von_mises_series_length(kappa));
  }
  return series.ratios;
}

inline double
VonMisesBicop::g(double theta, double concentration) const
{
  return std::exp(concentration * std::cos(theta)) /
         (tools_circular::two_pi() *
          boost::math::cyl_bessel_i(0, concentration));
}

//! \f$ \tilde G(\theta) = \bigl[\theta + 2 \sum_{j \ge 1} (I_j(\kappa) /
//! I_0(\kappa)) \sin(j\theta) / j\bigr] / 2\pi \f$.
inline double
VonMisesBicop::lifted_cdf(double theta, double concentration) const
{
  thread_local Series series;
  const auto& ratios = ratios_for(series, concentration);
  double sum = 0.0;
  for (size_t j = 1; j <= ratios.size(); ++j) {
    sum += ratios[j - 1] * std::sin(static_cast<double>(j) * theta) /
           static_cast<double>(j);
  }
  return (theta + 2.0 * sum) / tools_circular::two_pi();
}

//! the root lies within half a turn of \f$ 2\pi w \f$.
inline double
VonMisesBicop::lifted_cdf_inverse(double w, double concentration) const
{
  const double center = tools_circular::two_pi() * w;
  const double pi = boost::math::constants::pi<double>();
  return tools_circular::invert_increasing(
    [&](double t) { return lifted_cdf(t, concentration); },
    [&](double t) { return g(t, concentration); },
    w,
    center - pi,
    center + pi);
}

inline double
VonMisesBicop::concentration_from_resultant(double rbar) const
{
  return tools_circular::von_mises_a_inverse(rbar, parameters_upper_bounds_(0));
}
}
