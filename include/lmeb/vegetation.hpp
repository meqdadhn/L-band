// Vegetation modelling with the tau-omega model (Section 2.3).
#pragma once

#include "lmeb/common.hpp"

namespace lmeb {

struct VegetationParams {
  // Standing vegetation, Eqs. (22)-(24).
  double b_s1 = 0.06;   // b'_S  (slope of tau_S_NAD vs LAI)
  double b_s2 = 0.0;    // b''_S (intercept)
  double tt_h = 1.0;    // angular dependence of tau at H
  double tt_v = 1.0;    // angular dependence of tau at V
  // Single scattering albedo (independent of angle).
  double omega_h = 0.0;
  double omega_v = 0.0;
  // Composite temperature, Eq. (17).
  double bt = 1.7;
  // Litter, Eqs. (25)-(27).
  double c_l = 0.24;            // litter attenuation coefficient
  double a_l = 0.0;             // Mg_L = a_L * SM + b_L
  double b_l = 0.0;
  double litter_biomass = 0.0;  // BS_L, dry litter biomass (kg/m^2)
  // Interception flag threshold on polarization ratio at theta ~ 50 deg.
  double pr_interception = 0.02;
};

// tau_S_NAD = b'_S * LAI + b''_S, Eq. (22) (clipped at 0).
double standingNadirOpticalDepth(double lai, const VegetationParams& p);

// tau_SP(theta) = tau_S_NAD (sin^2(theta) tt_P + cos^2(theta)), Eqs. (23)-(24).
DualPol standingOpticalDepth(double tau_s_nad, double theta_rad, const VegetationParams& p);

// Mg_L = a_L * SM + b_L, clipped to [0, 0.8], Eq. (27).
double litterMoistureContent(double sm, const VegetationParams& p);

// LWC = Mg_L / (1 - Mg_L) * BS_L, Eq. (26).
double litterWaterContent(double mg_l, double litter_biomass);

// tau_L = c_L * LWC, Eq. (25).
double litterOpticalDepth(double sm, const VegetationParams& p);

// gamma_P = exp(-tau_P / cos(theta)), Eq. (15).
DualPol vegetationAttenuation(const DualPol& tau, double theta_rad);

// T_GC = A_t T_C + (1 - A_t) T_G, A_t = min(1, B_t (1 - exp(-tau_NAD))), Eqs. (16)-(18).
double compositeTemperature(double t_canopy, double t_ground, double tau_nad, double bt);

// Eq. (14): TB_P = (1-w)(1-g)(1+g r) T_C + (1-r) g T_G.
double tauOmegaBrightness(double omega, double gamma, double r_soil, double t_canopy,
                          double t_ground);

// PR = (TBV - TBH) / (TBV + TBH).
double polarizationRatio(const DualPol& tb);

// True when PR (at large incidence angles, theta >= 40 deg) is below the
// threshold PR_I, flagging significant rainfall interception.
bool interceptionFlag(const DualPol& tb, const VegetationParams& p);

}  // namespace lmeb
