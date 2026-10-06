#include "lmeb/vegetation.hpp"

#include <algorithm>

namespace lmeb {

double standingNadirOpticalDepth(double lai, const VegetationParams& p) {
  return std::max(0.0, p.b_s1 * lai + p.b_s2);
}

DualPol standingOpticalDepth(double tau_s_nad, double theta_rad, const VegetationParams& p) {
  const double s2 = std::pow(std::sin(theta_rad), 2);
  const double c2 = std::pow(std::cos(theta_rad), 2);
  return {tau_s_nad * (s2 * p.tt_h + c2), tau_s_nad * (s2 * p.tt_v + c2)};
}

double litterMoistureContent(double sm, const VegetationParams& p) {
  return std::clamp(p.a_l * sm + p.b_l, 0.0, 0.8);
}

double litterWaterContent(double mg_l, double litter_biomass) {
  return mg_l / (1.0 - mg_l) * litter_biomass;
}

double litterOpticalDepth(double sm, const VegetationParams& p) {
  if (p.litter_biomass <= 0.0) return 0.0;
  return p.c_l * litterWaterContent(litterMoistureContent(sm, p), p.litter_biomass);
}

DualPol vegetationAttenuation(const DualPol& tau, double theta_rad) {
  const double c = std::cos(theta_rad);
  return {std::exp(-tau.h / c), std::exp(-tau.v / c)};
}

double compositeTemperature(double t_canopy, double t_ground, double tau_nad, double bt) {
  const double at = std::min(1.0, bt * (1.0 - std::exp(-tau_nad)));
  return at * t_canopy + (1.0 - at) * t_ground;
}

double tauOmegaBrightness(double omega, double gamma, double r_soil, double t_canopy,
                          double t_ground) {
  return (1.0 - omega) * (1.0 - gamma) * (1.0 + gamma * r_soil) * t_canopy +
         (1.0 - r_soil) * gamma * t_ground;
}

double polarizationRatio(const DualPol& tb) { return (tb.v - tb.h) / (tb.v + tb.h); }

bool interceptionFlag(const DualPol& tb, const VegetationParams& p) {
  return polarizationRatio(tb) < p.pr_interception;
}

}  // namespace lmeb
