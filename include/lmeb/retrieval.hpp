// Multi-angular, dual-polarization retrieval (Section 3.2).
//
// Minimizes the cost function of Eq. (28):
//   CF = sum (TB_mes - TB*)^2 / sigma(TB)^2 + sum_i (P_i^ini - P_i*)^2 / sigma(P_i)^2
// with a Levenberg-Marquardt generalized least-squares scheme.
#pragma once

#include <string>
#include <vector>

#include "lmeb/model.hpp"

namespace lmeb {

enum class Param { SM, TauNad, TGC, HR, NRH, NRV, OmegaH, OmegaV, TtH, TtV };

std::string paramName(Param p);

struct Observation {
  double theta_deg;
  Pol pol;
  double tb;            // measured TB (K)
  double sigma = 2.0;   // sigma(TB) (K)
};

// A retrieved parameter with its a-priori value and constraint.
struct ParamPrior {
  Param param;
  double initial;
  double sigma;
};

struct RetrievalOptions {
  int max_iterations = 100;
  double tolerance = 1e-8;  // relative cost decrease to stop
  double initial_lambda = 1e-3;
};

struct RetrievalResult {
  std::vector<ParamPrior> priors;  // same order as the input priors
  std::vector<double> values;      // retrieved values P_i*
  std::vector<double> std_devs;    // posterior std. dev. from (J^T J)^-1
  double cost = 0.0;               // CF at the solution
  double rmse_tb = 0.0;            // RMSE between measured and simulated TB (K)
  int iterations = 0;
  bool converged = false;

  double value(Param p) const;  // throws if p was not retrieved
};

// Parameters not listed in `priors` are held fixed at their values in
// `config` / `state`.
RetrievalResult retrieve(const ModelConfig& config, const SurfaceState& state,
                         const std::vector<Observation>& observations,
                         const std::vector<ParamPrior>& priors,
                         const RetrievalOptions& options = {});

// Apply a parameter value to a configuration / state.
void setParam(Param p, double value, ModelConfig& config, SurfaceState& state);
double getParam(Param p, const ModelConfig& config, const SurfaceState& state);

// Simulated TB for each observation geometry.
std::vector<double> simulateObservations(const ModelConfig& config, const SurfaceState& state,
                                         const std::vector<Observation>& observations);

// The '3-P' configuration of Table 3: SM, tau_NAD and HR retrieved, T_GC fixed
// (sigma = 0.1 K) at the value given by Eqs. (16)-(17).
std::vector<ParamPrior> threeParamPriors(const ModelConfig& config, const SurfaceState& state,
                                         double tau_nad_initial, double hr_initial);

}  // namespace lmeb
