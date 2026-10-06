// L-MEB forward model for the nominal class (bare soil and low vegetation).
#pragma once

#include <limits>

#include "lmeb/common.hpp"
#include "lmeb/dielectric.hpp"
#include "lmeb/soil.hpp"
#include "lmeb/vegetation.hpp"

namespace lmeb {

// Static (calibrated / tabulated) parameters of one cover type.
struct ModelConfig {
  double freq_hz = kSmosFrequencyHz;
  SoilTexture texture;
  RoughnessParams roughness;
  EffectiveTempParams eff_temp;
  VegetationParams vegetation;
  // If true, the composite temperature T_GC (Eq. 16) replaces both T_C and
  // T_G in Eq. (14), as in the L-MEB retrievals. Otherwise T_C and T_G are used.
  bool use_composite_temperature = true;
};

// Time-varying surface state.
struct SurfaceState {
  double sm = 0.2;              // 0-3 cm volumetric soil moisture (m3/m3)
  double t_soil_surf = 295.0;   // K, ~0-5 cm
  double t_soil_depth = 290.0;  // K, ~50-100 cm
  double t_canopy = 293.0;      // K
  double tau_s_nad = 0.0;       // nadir optical depth of standing vegetation
  double tau_interception = 0.0;  // tau_IP (isotropic), Eq. (19)
  // Optional override of T_GC (used when T_GC is a retrieved parameter).
  double tgc = std::numeric_limits<double>::quiet_NaN();
};

struct ForwardResult {
  DualPol tb;            // surface-level brightness temperature (K)
  DualPol reflectivity;  // rough soil reflectivity r_GP
  DualPol tau;           // total optical depth tau_P, Eq. (19)
  DualPol gamma;         // vegetation attenuation gamma_P
  Complex eps_soil;
  double t_ground = 0.0;     // T_G
  double t_composite = 0.0;  // T_GC
  double tau_nad = 0.0;      // total nadir optical depth
};

// Simulate TB at incidence angle theta (radians).
ForwardResult simulate(const ModelConfig& config, const SurfaceState& state, double theta_rad);

}  // namespace lmeb
