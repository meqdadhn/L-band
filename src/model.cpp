#include "lmeb/model.hpp"

namespace lmeb {

ForwardResult simulate(const ModelConfig& cfg, const SurfaceState& st, double theta_rad) {
  ForwardResult out;

  // Soil.
  out.eps_soil = soilPermittivity(st.sm, st.t_soil_surf, cfg.texture, cfg.freq_hz);
  out.reflectivity =
      roughReflectivity(fresnelReflectivity(out.eps_soil, theta_rad), cfg.roughness, theta_rad);
  out.t_ground = effectiveSoilTemperature(st.sm, st.t_soil_surf, st.t_soil_depth, cfg.eff_temp);

  // Vegetation optical depth, Eq. (19): tau_P = tau_SP + tau_L + tau_IP.
  const double tau_isotropic = litterOpticalDepth(st.sm, cfg.vegetation) + st.tau_interception;
  const DualPol tau_s = standingOpticalDepth(st.tau_s_nad, theta_rad, cfg.vegetation);
  out.tau = {tau_s.h + tau_isotropic, tau_s.v + tau_isotropic};
  out.tau_nad = st.tau_s_nad + tau_isotropic;
  out.gamma = vegetationAttenuation(out.tau, theta_rad);

  out.t_composite = std::isnan(st.tgc)
                        ? compositeTemperature(st.t_canopy, out.t_ground, out.tau_nad,
                                               cfg.vegetation.bt)
                        : st.tgc;

  const double tc = cfg.use_composite_temperature ? out.t_composite : st.t_canopy;
  const double tg = cfg.use_composite_temperature ? out.t_composite : out.t_ground;

  out.tb.h = tauOmegaBrightness(cfg.vegetation.omega_h, out.gamma.h, out.reflectivity.h, tc, tg);
  out.tb.v = tauOmegaBrightness(cfg.vegetation.omega_v, out.gamma.v, out.reflectivity.v, tc, tg);
  return out;
}

}  // namespace lmeb
