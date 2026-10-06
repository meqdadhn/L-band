// '3-P' retrieval (SM, tau_NAD, HR) from multi-angular dual-polarization TB.
//
// Synthetic mode (default): simulates noisy observations from a known truth and
// retrieves it.
//   lmeb_retrieve --crop soybean --sm 0.22 --tau 0.15 --noise 2
//
// File mode: retrieves from a CSV file with lines "theta_deg,pol,tb[,sigma]"
// (pol is H or V, '#' starts a comment).
//   lmeb_retrieve --crop corn --obs obs.csv --tsurf 296 --tdepth 291 --tc 295
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

#include "cli.hpp"
#include "lmeb/lmeb.hpp"

namespace {

void usage() {
  std::cout <<
      "Usage: lmeb_retrieve [options]\n"
      "  --crop NAME        preset (default soybean)\n"
      "  --obs FILE         CSV observations theta,pol,tb[,sigma]; omit for synthetic mode\n"
      "  --tsurf K --tdepth K --tc K    temperatures (default 295/290/293)\n"
      "  --tau-init X       a-priori tau_NAD (default 0, Table 3)\n"
      "  --hr-init X        a-priori HR (default: preset value)\n"
      "  --sigma-tb K       TB uncertainty for synthetic data (default 2)\n"
      "  --tgc K            fixed composite temperature T_GC (default: Eqs. 16-17)\n"
      "  --sm-tg X          SM used for C_t when estimating T_GC in file mode (default 0.2)\n"
      " synthetic mode only:\n"
      "  --sm X --tau X --hr X   true SM, tau_NAD, HR (default 0.2, 0.15, preset)\n"
      "  --noise K          Gaussian noise std on TB (default 2)\n"
      "  --seed N           random seed (default 42)\n";
}

std::vector<lmeb::Observation> readObservations(const std::string& path) {
  std::ifstream in(path);
  if (!in) cli::Args::fail("cannot open " + path);
  std::vector<lmeb::Observation> obs;
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    if (const auto hash = line.find('#'); hash != std::string::npos) line.erase(hash);
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::stringstream ss(line);
    std::string theta, pol, tb, sigma;
    std::getline(ss, theta, ',');
    std::getline(ss, pol, ',');
    std::getline(ss, tb, ',');
    std::getline(ss, sigma, ',');
    try {
      lmeb::Observation o;
      o.theta_deg = std::stod(theta);
      pol.erase(0, pol.find_first_not_of(" \t"));
      if (pol.empty() || (pol[0] != 'H' && pol[0] != 'h' && pol[0] != 'V' && pol[0] != 'v'))
        throw std::invalid_argument("pol");
      o.pol = (pol[0] == 'H' || pol[0] == 'h') ? lmeb::Pol::H : lmeb::Pol::V;
      o.tb = std::stod(tb);
      if (!sigma.empty() && sigma.find_first_not_of(" \t\r") != std::string::npos)
        o.sigma = std::stod(sigma);
      obs.push_back(o);
    } catch (const std::exception&) {
      cli::Args::fail(path + ":" + std::to_string(lineno) + ": expected theta,pol,tb[,sigma]");
    }
  }
  return obs;
}

}  // namespace

int main(int argc, char** argv) {
  cli::Args args(argc, argv);
  if (args.has("help")) {
    usage();
    return 0;
  }

  lmeb::CropPreset preset;
  try {
    preset = lmeb::findCropPreset(args.str("crop", "soybean"));
  } catch (const std::exception& e) {
    cli::Args::fail(e.what());
  }
  const lmeb::ModelConfig cfg = lmeb::makeConfig(preset);

  lmeb::SurfaceState st;
  st.t_soil_surf = args.num("tsurf", 295.0);
  st.t_soil_depth = args.num("tdepth", 290.0);
  st.t_canopy = args.num("tc", 293.0);
  const double sigma_tb = args.num("sigma-tb", 2.0);

  std::vector<lmeb::Observation> obs;
  const bool synthetic = !args.has("obs");
  lmeb::ModelConfig truth_cfg = cfg;
  lmeb::SurfaceState truth = st;
  if (synthetic) {
    truth.sm = args.num("sm", 0.2);
    truth.tau_s_nad = args.num("tau", 0.15);
    truth_cfg.roughness.hr = args.num("hr", cfg.roughness.hr);
    std::mt19937 rng(static_cast<unsigned>(args.num("seed", 42)));
    std::normal_distribution<double> noise(0.0, args.num("noise", 2.0));
    for (double th = 10.0; th <= 55.0; th += 5.0) {  // theta <= 55 deg as in the paper
      const lmeb::DualPol tb = lmeb::simulate(truth_cfg, truth, lmeb::deg2rad(th)).tb;
      obs.push_back({th, lmeb::Pol::H, tb.h + noise(rng), sigma_tb});
      obs.push_back({th, lmeb::Pol::V, tb.v + noise(rng), sigma_tb});
    }
  } else {
    obs = readObservations(args.str("obs", ""));
    if (obs.empty()) cli::Args::fail("no observations read");
  }

  // T_GC is held fixed (sigma = 0.1 K, Table 3). It is estimated from Eqs. (11)-(12)
  // and (16)-(17), which need an SM estimate for C_t; in synthetic mode the truth is used.
  lmeb::SurfaceState prior_state = st;
  prior_state.sm = synthetic ? truth.sm : args.num("sm-tg", 0.2);
  auto priors = lmeb::threeParamPriors(cfg, prior_state, args.num("tau-init", 0.0),
                                       args.num("hr-init", cfg.roughness.hr));
  if (args.has("tgc"))
    for (auto& p : priors)
      if (p.param == lmeb::Param::TGC) p.initial = args.num("tgc", p.initial);
  const lmeb::RetrievalResult res = lmeb::retrieve(cfg, st, obs, priors);

  std::printf("# L-MEB '3-P' retrieval: %s, %zu observations\n", preset.name.c_str(), obs.size());
  std::printf("# converged=%s iterations=%d cost=%.4f RMSE(TB)=%.3f K\n",
              res.converged ? "yes" : "no", res.iterations, res.cost, res.rmse_tb);
  std::printf("%-10s %10s %10s %10s", "param", "initial", "retrieved", "std");
  if (synthetic) std::printf(" %10s", "truth");
  std::printf("\n");
  for (std::size_t i = 0; i < res.priors.size(); ++i) {
    const lmeb::Param p = res.priors[i].param;
    std::printf("%-10s %10.4f %10.4f %10.4f", lmeb::paramName(p).c_str(), res.priors[i].initial,
                res.values[i], res.std_devs[i]);
    if (synthetic) {
      double t = lmeb::getParam(p, truth_cfg, truth);
      if (p == lmeb::Param::TGC) t = lmeb::simulate(truth_cfg, truth, 0.0).t_composite;
      std::printf(" %10.4f", t);
    }
    std::printf("\n");
  }
  return 0;
}
