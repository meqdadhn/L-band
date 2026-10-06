// Soil emission: reflectivity (Section 2.2.1) and effective temperature (Section 2.2.2).
#pragma once

#include <vector>

#include "lmeb/common.hpp"

namespace lmeb {

// Fresnel reflectivity r*_GP = |R_P(eps_G, theta)|^2 of a flat surface, Eq. (3).
DualPol fresnelReflectivity(Complex eps, double theta_rad);

// Coherent reflectivity of a Gaussian rough surface, Eq. (5).
// sd_m: standard deviation of surface height (m).
DualPol coherentReflectivity(const DualPol& fresnel, double sd_m, double theta_rad,
                             double freq_hz = kSmosFrequencyHz);

// Semi-empirical roughness parameters of Eq. (7).
struct RoughnessParams {
  double hr = 0.1;   // H_R
  double nr_h = 0.0; // N_RH
  double nr_v = 0.0; // N_RV
  double qr = 0.0;   // Q_R (polarization mixing; ~0 at L-band)
};

// Rough-surface reflectivity r_GP, Eq. (7) (reduces to Eq. (8) for Q_R = 0):
//   r_GP = [(1-Q_R) r*_GP + Q_R r*_GQ] exp(-H_R cos^{N_RP}(theta))
DualPol roughReflectivity(const DualPol& fresnel, const RoughnessParams& rough, double theta_rad);

// Parameters of the effective-temperature parameterization, Eq. (12).
struct EffectiveTempParams {
  double w0 = 0.3;   // m3/m3
  double bw0 = 0.3;  // default values used in L-MEB
};

// T_G = T_depth + C_t (T_surf - T_depth), Eq. (11), with constant C_t
// (Choudhury et al. 1982 found C_t = 0.246 at L-band).
double effectiveSoilTemperature(double t_surf, double t_depth, double ct);

// C_t = (SM / w0)^bw0, Eq. (12), clipped to [0, 1].
double effectiveTempCoefficient(double sm, const EffectiveTempParams& p);

// T_G from Eqs. (11)-(12) (the L-MEB default).
double effectiveSoilTemperature(double sm, double t_surf, double t_depth,
                                const EffectiveTempParams& p);

// One layer of a soil profile, from the surface downwards.
struct SoilLayer {
  double thickness_m;  // layer thickness; the last layer is treated as semi-infinite
  double temp_k;
  Complex eps;
};

// Radiative-transfer effective temperature, Eqs. (9)-(10), evaluated exactly
// for a piecewise-constant profile.
double effectiveSoilTemperature(const std::vector<SoilLayer>& profile,
                                double freq_hz = kSmosFrequencyHz);

}  // namespace lmeb
