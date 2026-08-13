// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Dumps a JSON file of numerical outputs across the surfaces touched by the
// performance work, for before/after comparison via
// scripts/compare_parity.py. Fully deterministic (fixed seeds everywhere).
//
// Usage: parity_dump [output.json]  (default: bench_results/parity.json)
// Keep generated dumps under bench_results/ so they stay gitignored.

#include "helpers.hpp"
#include <fstream>
#include <iostream>
#include <vinecopulib.hpp>
#include <vinecopulib/misc/nlohmann_json.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

using namespace vinecopulib;
using json = nlohmann::json;

namespace {

std::vector<double>
to_vec(const Eigen::MatrixXd& m)
{
  return std::vector<double>(m.data(), m.data() + m.size());
}

//! 21x21 lattice on (0.025, ..., 0.975)^2.
Eigen::MatrixXd
lattice()
{
  Eigen::VectorXd g = Eigen::VectorXd::LinSpaced(21, 0.025, 0.975);
  Eigen::MatrixXd u(21 * 21, 2);
  size_t k = 0;
  for (long i = 0; i < 21; ++i) {
    for (long j = 0; j < 21; ++j) {
      u(k, 0) = g(i);
      u(k, 1) = g(j);
      ++k;
    }
  }
  return u;
}

json
dump_bicop_eval()
{
  json out;
  const auto u = lattice();
  const std::vector<std::pair<BicopFamily, int>> cases = {
    { BicopFamily::gaussian, 0 }, { BicopFamily::student, 0 },
    { BicopFamily::clayton, 0 },  { BicopFamily::clayton, 90 },
    { BicopFamily::gumbel, 0 },   { BicopFamily::frank, 0 },
    { BicopFamily::joe, 0 },      { BicopFamily::bb1, 0 },
    { BicopFamily::bb6, 0 },      { BicopFamily::bb7, 0 },
    { BicopFamily::bb8, 0 },      { BicopFamily::tawn, 0 },
  };
  for (const auto& c : cases) {
    Bicop bc(c.first, c.second, bench::family_parameters(c.first));
    const std::string key =
      get_family_name(c.first) + "/rot" + std::to_string(c.second);
    out[key]["pdf"] = to_vec(bc.pdf(u));
    out[key]["cdf"] = to_vec(bc.cdf(u));
    out[key]["hfunc1"] = to_vec(bc.hfunc1(u));
    out[key]["hfunc2"] = to_vec(bc.hfunc2(u));
    out[key]["hinv1"] = to_vec(bc.hinv1(u));
    out[key]["hinv2"] = to_vec(bc.hinv2(u));
    out[key]["loglik"] = bc.loglik(u);
    out[key]["tau"] = bc.parameters_to_tau(bc.get_parameters());
  }
  // discrete dispatch
  for (auto family : { BicopFamily::gaussian, BicopFamily::clayton }) {
    Bicop bc(family, 0, bench::family_parameters(family));
    bc.set_var_types({ "d", "c" });
    const auto u_disc = bench::discretize_first(u);
    const std::string key = get_family_name(family) + "/disc";
    out[key]["pdf"] = to_vec(bc.pdf(u_disc));
    out[key]["hfunc1"] = to_vec(bc.hfunc1(u_disc));
    out[key]["hfunc2"] = to_vec(bc.hfunc2(u_disc));
    out[key]["hinv1"] = to_vec(bc.hinv1(u_disc));
  }
  return out;
}

json
dump_bicop_fits()
{
  json out;
  const std::vector<BicopFamily> families = {
    BicopFamily::gaussian, BicopFamily::student, BicopFamily::clayton,
    BicopFamily::gumbel,   BicopFamily::frank,   BicopFamily::joe,
    BicopFamily::bb1,      BicopFamily::bb6,     BicopFamily::bb7,
    BicopFamily::bb8,      BicopFamily::tawn
  };
  for (auto family : families) {
    const auto data = bench::sim_data(family, 0, 500);
    const std::string fam_name = get_family_name(family);
    {
      FitControlsBicop controls({ family });
      controls.set_parametric_method("mle");
      Bicop bc(family);
      bc.fit(data, controls);
      out["mle/" + fam_name]["par"] = to_vec(bc.get_parameters());
      out["mle/" + fam_name]["loglik"] = bc.get_loglik();
    }
    const auto& itau = bicop_families::itau;
    if (std::find(itau.begin(), itau.end(), family) != itau.end()) {
      FitControlsBicop controls({ family });
      controls.set_parametric_method("itau");
      Bicop bc(family);
      bc.fit(data, controls);
      out["itau/" + fam_name]["par"] = to_vec(bc.get_parameters());
      out["itau/" + fam_name]["loglik"] = bc.get_loglik();
    }
  }
  return out;
}

json
dump_tll()
{
  json out;
  const auto data = bench::sim_data(BicopFamily::gaussian, 0, 500);
  const auto u = lattice();
  for (const auto& method : { std::string("constant"),
                              std::string("linear"),
                              std::string("quadratic") }) {
    FitControlsBicop controls({ BicopFamily::tll });
    controls.set_nonparametric_method(method);
    Bicop bc(BicopFamily::tll);
    bc.fit(data, controls);
    json& j = out[method];
    const Eigen::MatrixXd grid = bc.get_parameters();
    // probe every 37th grid value (900 -> 25 probes) + full shape
    std::vector<double> probes;
    for (long k = 0; k < grid.size(); k += 37) {
      probes.push_back(grid(k));
    }
    j["grid_probes"] = probes;
    j["npars"] = bc.get_npars();
    j["loglik"] = bc.get_loglik();
    j["pdf"] = to_vec(bc.pdf(u));
    j["hfunc1"] = to_vec(bc.hfunc1(u));
    j["hfunc2"] = to_vec(bc.hfunc2(u));
    j["hinv1"] = to_vec(bc.hinv1(u));
    j["cdf"] = to_vec(bc.cdf(u.topRows(21)));
  }
  // discrete tll
  {
    const auto data_disc = bench::discretize_first(data);
    FitControlsBicop controls({ BicopFamily::tll });
    Bicop bc(BicopFamily::tll);
    bc.set_var_types({ "d", "c" });
    bc.fit(data_disc, controls);
    out["disc_constant"]["npars"] = bc.get_npars();
    out["disc_constant"]["loglik"] = bc.get_loglik();
    out["disc_constant"]["pdf"] = to_vec(bc.pdf(bench::discretize_first(u)));
  }
  return out;
}

json
dump_tools_stats()
{
  json out;
  const auto u = lattice();
  const Eigen::MatrixXd z = tools_stats::qnorm(u);
  for (double rho : { 0.1, 0.5, 0.9 }) {
    out["pbvnorm/rho=" + std::to_string(rho).substr(0, 3)] =
      to_vec(tools_stats::pbvnorm(z, rho));
  }
  const Eigen::MatrixXd zt4 = tools_stats::qt(u, 4.0);
  const Eigen::MatrixXd zt5 = tools_stats::qt(u, 5.0);
  out["pbvt/nu=4"] = to_vec(tools_stats::pbvt(zt4, 4, 0.5));
  out["pbvt/nu=5"] = to_vec(tools_stats::pbvt(zt5, 5, 0.5));

  const auto x = tools_stats::simulate_uniform(50, 3, false, { 5 });
  out["pseudo_obs/average"] = to_vec(tools_stats::to_pseudo_obs(x));
  Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(50, 0.5, 1.5);
  out["pseudo_obs/weighted"] =
    to_vec(tools_stats::to_pseudo_obs(x, "average", w));
  out["pseudo_obs/random"] =
    to_vec(tools_stats::to_pseudo_obs(x, "random", Eigen::VectorXd(), { 17 }));

  out["ghalton"] = to_vec(tools_stats::ghalton(8, 5, { 5 }));
  out["sobol"] = to_vec(tools_stats::sobol(8, 5, { 5 }));
  out["simulate_uniform"] =
    to_vec(tools_stats::simulate_uniform(8, 3, false, { 5 }));
  out["simulate_uniform_qrng"] =
    to_vec(tools_stats::simulate_uniform(8, 3, true, { 5 }));

  const auto xm = bench::sim_data(BicopFamily::gaussian, 0, 500);
  out["pairwise_mcor"] = tools_stats::pairwise_mcor(xm);
  return out;
}

json
dump_tau_sweeps()
{
  json out;
  auto sweep = [](BicopFamily family,
                  const std::vector<Eigen::VectorXd>& pars) {
    Bicop bc(family);
    std::vector<double> taus;
    for (const auto& p : pars) {
      taus.push_back(bc.parameters_to_tau(p));
    }
    return taus;
  };
  std::vector<Eigen::VectorXd> bb_pars;
  for (double th : { 1.5, 2.0, 3.0 }) {
    for (double de : { 1.5, 2.5 }) {
      Eigen::VectorXd p(2);
      p << th, de;
      bb_pars.push_back(p);
    }
  }
  out["bb6"] = sweep(BicopFamily::bb6, bb_pars);
  out["bb7"] = sweep(BicopFamily::bb7, bb_pars);
  {
    std::vector<Eigen::VectorXd> pars;
    for (double th : { 2.0, 4.0 }) {
      for (double de : { 0.4, 0.8 }) {
        Eigen::VectorXd p(2);
        p << th, de;
        pars.push_back(p);
      }
    }
    out["bb8"] = sweep(BicopFamily::bb8, pars);
  }
  {
    std::vector<Eigen::VectorXd> pars;
    for (double psi1 : { 0.3, 0.8 }) {
      for (double th : { 1.5, 3.0 }) {
        Eigen::VectorXd p(3);
        p << psi1, 0.5, th;
        pars.push_back(p);
      }
    }
    out["tawn"] = sweep(BicopFamily::tawn, pars);
  }
  // Near-boundary corners, kept in separate keys so the interior sweeps above
  // stay index-stable. The integrands are worst-conditioned here, so this is
  // where a quadrature or reformulation change shows up first.
  auto corners = [](const std::vector<std::vector<double>>& rows) {
    std::vector<Eigen::VectorXd> pars;
    for (const auto& row : rows) {
      Eigen::VectorXd p(row.size());
      for (size_t i = 0; i < row.size(); ++i) {
        p(static_cast<Eigen::Index>(i)) = row[i];
      }
      pars.push_back(p);
    }
    return pars;
  };
  out["bb6_edge"] = sweep(
    BicopFamily::bb6,
    corners({ { 1.0, 1.0 }, { 1.001, 8.0 }, { 6.0, 1.0 }, { 6.0, 8.0 } }));
  out["bb7_edge"] = sweep(
    BicopFamily::bb7,
    corners({ { 1.0, 0.01 }, { 1.001, 25.0 }, { 6.0, 0.01 }, { 6.0, 25.0 } }));
  out["bb8_edge"] = sweep(
    BicopFamily::bb8,
    corners({ { 1.0, 1e-4 }, { 1.001, 1.0 }, { 8.0, 1e-4 }, { 8.0, 1.0 } }));
  out["tawn_edge"] = sweep(BicopFamily::tawn,
                           corners({ { 0.001, 0.5, 60.0 },
                                     { 0.5, 0.5, 50.0 },
                                     { 0.5, 0.5, 60.0 },
                                     { 0.999, 0.999, 60.0 },
                                     { 1.0, 1.0, 60.0 } }));
  return out;
}

json
dump_vinecop()
{
  json out;
  auto vc = bench::make_gaussian_vine(5);
  const auto u = vc.simulate(100, false, 1, { 5 });
  out["simulate"] = to_vec(u.topRows(3));
  out["pdf"] = to_vec(vc.pdf(u));
  out["loglik"] = vc.loglik(u);
  out["rosenblatt"] = to_vec(vc.rosenblatt(u).topRows(5));
  out["inverse_rosenblatt"] = to_vec(vc.inverse_rosenblatt(u).topRows(5));
  out["cdf"] = to_vec(vc.cdf(u.topRows(5), 1000, 1, { 5 }));
  out["scores_stepwise"] = to_vec(vc.scores(u, true, 1).colwise().sum());
  out["scores_full"] = to_vec(vc.scores(u, false, 1).colwise().sum());
  out["hessian_avg_stepwise"] = to_vec(vc.hessian(u, true, 1));
  out["hessian_avg_joint"] = to_vec(vc.hessian(u, false, 1));
  out["scores_cov"] = to_vec(vc.scores_cov(u, true, 1));

  // selection determinism (itau, fixed seed data)
  FitControlsVinecop controls(bicop_families::itau, "itau");
  Vinecop fitted(u, RVineStructure(), {}, controls);
  out["select/matrix"] = [&] {
    const auto m = fitted.get_matrix();
    std::vector<double> v;
    for (long i = 0; i < m.rows(); ++i)
      for (long j = 0; j < m.cols(); ++j)
        v.push_back(static_cast<double>(m(i, j)));
    return v;
  }();
  out["select/loglik"] = fitted.get_loglik();
  json pars = json::array();
  for (const auto& tree : fitted.get_all_parameters()) {
    for (const auto& p : tree) {
      pars.push_back(to_vec(p));
    }
  }
  out["select/parameters"] = pars;
  return out;
}

} // namespace

int
main(int argc, char** argv)
{
  const std::string path = (argc > 1) ? argv[1] : "bench_results/parity.json";
  json j;
  std::cerr << "bicop_eval..." << std::endl;
  j["bicop_eval"] = dump_bicop_eval();
  std::cerr << "bicop_fit..." << std::endl;
  j["bicop_fit"] = dump_bicop_fits();
  std::cerr << "tll..." << std::endl;
  j["tll"] = dump_tll();
  std::cerr << "tools_stats..." << std::endl;
  j["tools_stats"] = dump_tools_stats();
  std::cerr << "tau..." << std::endl;
  j["tau"] = dump_tau_sweeps();
  std::cerr << "vinecop..." << std::endl;
  j["vinecop"] = dump_vinecop();

  std::ofstream file(path);
  if (!file.is_open()) {
    std::cerr << "could not open '" << path
              << "' for writing (does its directory exist?)" << std::endl;
    return 1;
  }
  file << j.dump(1);
  file.close();
  std::cout << "wrote " << path << std::endl;
  return 0;
}
