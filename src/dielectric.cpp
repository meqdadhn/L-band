#include "lmeb/dielectric.hpp"

#include <algorithm>

namespace lmeb {

namespace {
constexpr double kSoilSpecificDensity = 2.66;  // g/cm^3, density of soil solids
constexpr double kDobsonAlpha = 0.65;
constexpr double kWaterHighFreqPermittivity = 4.9;
}  // namespace

Complex freeWaterPermittivity(double freq_hz, double temp_c, double salinity_ppt) {
  const double t = temp_c;
  const double s = salinity_ppt;

  // Pure water static permittivity and relaxation time (Ulaby et al., 1986).
  const double eps_w0 = 88.045 - 0.4147 * t + 6.295e-4 * t * t + 1.075e-5 * t * t * t;
  const double two_pi_tau =
      1.1109e-10 - 3.824e-12 * t + 6.938e-14 * t * t - 5.096e-16 * t * t * t;

  // Salinity corrections (Klein & Swift, 1977).
  const double a = 1.0 + 1.613e-5 * t * s - 3.656e-3 * s + 3.210e-5 * s * s - 4.232e-7 * s * s * s;
  const double b = 1.0 + 2.282e-5 * t * s - 7.638e-4 * s - 7.760e-6 * s * s + 1.105e-8 * s * s * s;

  const double eps_static = eps_w0 * a;
  const double ft = freq_hz * two_pi_tau * b;
  const double denom = 1.0 + ft * ft;

  const double re = kWaterHighFreqPermittivity + (eps_static - kWaterHighFreqPermittivity) / denom;
  const double im = ft * (eps_static - kWaterHighFreqPermittivity) / denom;
  return {re, im};
}

Complex dobsonPermittivity(double sm, double temp_c, const SoilTexture& tx, double freq_hz) {
  sm = std::max(sm, 0.0);
  const double rho_b = tx.bulk_density;
  const double rho_s = kSoilSpecificDensity;

  const double eps_solid = std::pow(1.01 + 0.44 * rho_s, 2) - 0.062;
  const double beta_re = 1.2748 - 0.519 * tx.sand - 0.152 * tx.clay;
  const double beta_im = 1.33797 - 0.603 * tx.sand - 0.166 * tx.clay;
  // Effective conductivity (S/m), Dobson et al. (1985) / Peplinski et al. (1995), 1.4-18 GHz.
  const double sigma_eff = -1.645 + 1.939 * rho_b - 2.25622 * tx.sand + 1.594 * tx.clay;

  Complex eps_fw = freeWaterPermittivity(freq_hz, temp_c, tx.salinity);
  double eps_fw_im = eps_fw.imag();
  if (sm > 0.0) {
    eps_fw_im += sigma_eff / (2.0 * kPi * kVacuumPermittivity * freq_hz) * (rho_s - rho_b) /
                 (rho_s * sm);
  }

  const double re_alpha = 1.0 + rho_b / rho_s * (std::pow(eps_solid, kDobsonAlpha) - 1.0) +
                          std::pow(sm, beta_re) * std::pow(eps_fw.real(), kDobsonAlpha) - sm;
  const double im_alpha = std::pow(sm, beta_im) * std::pow(std::max(eps_fw_im, 0.0), kDobsonAlpha);

  return {std::pow(re_alpha, 1.0 / kDobsonAlpha), std::pow(im_alpha, 1.0 / kDobsonAlpha)};
}

Complex soilPermittivity(double sm, double temp_k, const SoilTexture& texture, double freq_hz) {
  if (temp_k < kZeroCelsius) return kFrozenSoilPermittivity;
  return dobsonPermittivity(sm, temp_k - kZeroCelsius, texture, freq_hz);
}

}  // namespace lmeb
