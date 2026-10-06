#include "lmeb/retrieval.hpp"

#include <algorithm>
#include <stdexcept>

namespace lmeb {

namespace {

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct Bounds {
  double lo, hi;
};

Bounds paramBounds(Param p) {
  switch (p) {
    case Param::SM: return {0.0, 0.7};
    case Param::TauNad: return {0.0, 5.0};
    case Param::TGC: return {200.0, 350.0};
    case Param::HR: return {0.0, 5.0};
    case Param::NRH:
    case Param::NRV: return {-3.0, 3.0};
    case Param::OmegaH:
    case Param::OmegaV: return {0.0, 0.5};
    case Param::TtH:
    case Param::TtV: return {0.0, 20.0};
  }
  return {-1e300, 1e300};
}

// Solves A x = b by Gaussian elimination with partial pivoting. Returns false if singular.
bool solve(Matrix a, Vector b, Vector& x) {
  const std::size_t n = b.size();
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    for (std::size_t r = col + 1; r < n; ++r)
      if (std::abs(a[r][col]) > std::abs(a[pivot][col])) pivot = r;
    if (std::abs(a[pivot][col]) < 1e-300) return false;
    std::swap(a[col], a[pivot]);
    std::swap(b[col], b[pivot]);
    for (std::size_t r = col + 1; r < n; ++r) {
      const double f = a[r][col] / a[col][col];
      for (std::size_t c = col; c < n; ++c) a[r][c] -= f * a[col][c];
      b[r] -= f * b[col];
    }
  }
  x.assign(n, 0.0);
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t c = i + 1; c < n; ++c) s -= a[i][c] * x[c];
    x[i] = s / a[i][i];
  }
  return true;
}

class Problem {
 public:
  Problem(const ModelConfig& cfg, const SurfaceState& st, const std::vector<Observation>& obs,
          const std::vector<ParamPrior>& priors)
      : cfg_(cfg), st_(st), obs_(obs), priors_(priors) {}

  std::size_t numResiduals() const { return obs_.size() + priors_.size(); }

  Vector residuals(const Vector& x) const {
    ModelConfig cfg = cfg_;
    SurfaceState st = st_;
    for (std::size_t i = 0; i < x.size(); ++i) setParam(priors_[i].param, x[i], cfg, st);
    const Vector sim = simulateObservations(cfg, st, obs_);
    Vector r(numResiduals());
    for (std::size_t k = 0; k < obs_.size(); ++k) r[k] = (obs_[k].tb - sim[k]) / obs_[k].sigma;
    for (std::size_t i = 0; i < priors_.size(); ++i)
      r[obs_.size() + i] = (priors_[i].initial - x[i]) / priors_[i].sigma;
    return r;
  }

  Matrix jacobian(const Vector& x, const Vector& r0) const {
    Matrix j(numResiduals(), Vector(x.size()));
    for (std::size_t i = 0; i < x.size(); ++i) {
      const Bounds b = paramBounds(priors_[i].param);
      double h = 1e-6 * std::max(1.0, std::abs(x[i]));
      if (x[i] + h > b.hi) h = -h;  // step inward at the upper bound
      Vector xp = x;
      xp[i] += h;
      const Vector rp = residuals(xp);
      for (std::size_t k = 0; k < r0.size(); ++k) j[k][i] = (rp[k] - r0[k]) / h;
    }
    return j;
  }

 private:
  const ModelConfig& cfg_;
  const SurfaceState& st_;
  const std::vector<Observation>& obs_;
  const std::vector<ParamPrior>& priors_;
};

double sumSquares(const Vector& r) {
  double s = 0.0;
  for (double v : r) s += v * v;
  return s;
}

void normalEquations(const Matrix& j, const Vector& r, Matrix& jtj, Vector& jtr) {
  const std::size_t n = j.empty() ? 0 : j[0].size();
  jtj.assign(n, Vector(n, 0.0));
  jtr.assign(n, 0.0);
  for (std::size_t k = 0; k < j.size(); ++k)
    for (std::size_t a = 0; a < n; ++a) {
      jtr[a] += j[k][a] * r[k];
      for (std::size_t b = 0; b < n; ++b) jtj[a][b] += j[k][a] * j[k][b];
    }
}

}  // namespace

std::string paramName(Param p) {
  switch (p) {
    case Param::SM: return "SM";
    case Param::TauNad: return "tau_NAD";
    case Param::TGC: return "T_GC";
    case Param::HR: return "HR";
    case Param::NRH: return "NRH";
    case Param::NRV: return "NRV";
    case Param::OmegaH: return "omega_H";
    case Param::OmegaV: return "omega_V";
    case Param::TtH: return "tt_H";
    case Param::TtV: return "tt_V";
  }
  return "?";
}

void setParam(Param p, double value, ModelConfig& cfg, SurfaceState& st) {
  switch (p) {
    case Param::SM: st.sm = value; break;
    case Param::TauNad: st.tau_s_nad = value; break;
    case Param::TGC: st.tgc = value; break;
    case Param::HR: cfg.roughness.hr = value; break;
    case Param::NRH: cfg.roughness.nr_h = value; break;
    case Param::NRV: cfg.roughness.nr_v = value; break;
    case Param::OmegaH: cfg.vegetation.omega_h = value; break;
    case Param::OmegaV: cfg.vegetation.omega_v = value; break;
    case Param::TtH: cfg.vegetation.tt_h = value; break;
    case Param::TtV: cfg.vegetation.tt_v = value; break;
  }
}

double getParam(Param p, const ModelConfig& cfg, const SurfaceState& st) {
  switch (p) {
    case Param::SM: return st.sm;
    case Param::TauNad: return st.tau_s_nad;
    case Param::TGC: return st.tgc;
    case Param::HR: return cfg.roughness.hr;
    case Param::NRH: return cfg.roughness.nr_h;
    case Param::NRV: return cfg.roughness.nr_v;
    case Param::OmegaH: return cfg.vegetation.omega_h;
    case Param::OmegaV: return cfg.vegetation.omega_v;
    case Param::TtH: return cfg.vegetation.tt_h;
    case Param::TtV: return cfg.vegetation.tt_v;
  }
  return 0.0;
}

double RetrievalResult::value(Param p) const {
  for (std::size_t i = 0; i < priors.size(); ++i)
    if (priors[i].param == p) return values[i];
  throw std::invalid_argument("parameter not retrieved: " + paramName(p));
}

std::vector<double> simulateObservations(const ModelConfig& cfg, const SurfaceState& st,
                                         const std::vector<Observation>& obs) {
  std::vector<double> tb(obs.size());
  for (std::size_t k = 0; k < obs.size(); ++k)
    tb[k] = simulate(cfg, st, deg2rad(obs[k].theta_deg)).tb[obs[k].pol];
  return tb;
}

std::vector<ParamPrior> threeParamPriors(const ModelConfig& cfg, const SurfaceState& st,
                                         double tau_nad_initial, double hr_initial) {
  SurfaceState s = st;
  s.tau_s_nad = tau_nad_initial;
  s.tgc = std::numeric_limits<double>::quiet_NaN();
  const double tgc = simulate(cfg, s, 0.0).t_composite;
  return {
      {Param::SM, 0.05, 0.3},
      {Param::TauNad, tau_nad_initial, 0.05},
      {Param::TGC, tgc, 0.1},
      {Param::HR, hr_initial, 0.1},
  };
}

RetrievalResult retrieve(const ModelConfig& cfg, const SurfaceState& st,
                         const std::vector<Observation>& obs,
                         const std::vector<ParamPrior>& priors, const RetrievalOptions& opt) {
  if (obs.empty()) throw std::invalid_argument("retrieve: no observations");
  for (const auto& p : priors)
    if (p.sigma <= 0.0) throw std::invalid_argument("retrieve: sigma must be > 0 for " + paramName(p.param));

  const Problem problem(cfg, st, obs, priors);
  const std::size_t n = priors.size();

  Vector x(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Bounds b = paramBounds(priors[i].param);
    x[i] = std::clamp(priors[i].initial, b.lo, b.hi);
  }

  Vector r = problem.residuals(x);
  double cost = sumSquares(r);
  double lambda = opt.initial_lambda;

  RetrievalResult res;
  res.priors = priors;

  Matrix jtj;
  Vector jtr;
  int it = 0;
  for (; it < opt.max_iterations; ++it) {
    normalEquations(problem.jacobian(x, r), r, jtj, jtr);

    bool improved = false;
    while (lambda < 1e12) {
      Matrix a = jtj;
      for (std::size_t i = 0; i < n; ++i) a[i][i] += lambda * std::max(jtj[i][i], 1e-12);
      Vector neg_jtr(n), delta;
      for (std::size_t i = 0; i < n; ++i) neg_jtr[i] = -jtr[i];
      if (!solve(a, neg_jtr, delta)) {
        lambda *= 10.0;
        continue;
      }
      Vector x_new(n);
      for (std::size_t i = 0; i < n; ++i) {
        const Bounds b = paramBounds(priors[i].param);
        x_new[i] = std::clamp(x[i] + delta[i], b.lo, b.hi);
      }
      const Vector r_new = problem.residuals(x_new);
      const double cost_new = sumSquares(r_new);
      if (cost_new < cost) {
        const double rel = (cost - cost_new) / std::max(cost, 1e-300);
        x = x_new;
        r = r_new;
        cost = cost_new;
        lambda = std::max(lambda / 10.0, 1e-12);
        improved = true;
        if (rel < opt.tolerance) res.converged = true;
        break;
      }
      lambda *= 10.0;
    }
    if (!improved) {
      res.converged = true;  // no descent direction left: at a (local) minimum
      break;
    }
    if (res.converged) break;
  }

  // Posterior covariance (J^T J)^-1 at the solution.
  normalEquations(problem.jacobian(x, r), r, jtj, jtr);
  res.std_devs.assign(n, std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < n; ++i) {
    Vector e(n, 0.0), col;
    e[i] = 1.0;
    if (solve(jtj, e, col) && col[i] >= 0.0) res.std_devs[i] = std::sqrt(col[i]);
  }

  double sq = 0.0;
  for (std::size_t k = 0; k < obs.size(); ++k) sq += std::pow(r[k] * obs[k].sigma, 2);
  res.rmse_tb = std::sqrt(sq / obs.size());
  res.values = x;
  res.cost = cost;
  res.iterations = it + 1;
  return res;
}

}  // namespace lmeb
