#include "lmeb/soil.hpp"

#include <algorithm>
#include <stdexcept>

namespace lmeb {

DualPol fresnelReflectivity(Complex eps, double theta_rad) {
  const double c = std::cos(theta_rad);
  const double s = std::sin(theta_rad);
  const Complex root = std::sqrt(eps - s * s);
  const Complex rh = (c - root) / (c + root);
  const Complex rv = (eps * c - root) / (eps * c + root);
  return {std::norm(rh), std::norm(rv)};
}

DualPol coherentReflectivity(const DualPol& fresnel, double sd_m, double theta_rad,
                             double freq_hz) {
  const double k = 4.0 * kPi * sd_m * std::cos(theta_rad) / wavelength(freq_hz);
  const double f = std::exp(-k * k);
  return {fresnel.h * f, fresnel.v * f};
}

DualPol roughReflectivity(const DualPol& fresnel, const RoughnessParams& r, double theta_rad) {
  const double c = std::cos(theta_rad);
  const double mixed_h = (1.0 - r.qr) * fresnel.h + r.qr * fresnel.v;
  const double mixed_v = (1.0 - r.qr) * fresnel.v + r.qr * fresnel.h;
  return {mixed_h * std::exp(-r.hr * std::pow(c, r.nr_h)),
          mixed_v * std::exp(-r.hr * std::pow(c, r.nr_v))};
}

double effectiveSoilTemperature(double t_surf, double t_depth, double ct) {
  return t_depth + ct * (t_surf - t_depth);
}

double effectiveTempCoefficient(double sm, const EffectiveTempParams& p) {
  if (sm <= 0.0) return 0.0;
  return std::clamp(std::pow(sm / p.w0, p.bw0), 0.0, 1.0);
}

double effectiveSoilTemperature(double sm, double t_surf, double t_depth,
                                const EffectiveTempParams& p) {
  return effectiveSoilTemperature(t_surf, t_depth, effectiveTempCoefficient(sm, p));
}

double effectiveSoilTemperature(const std::vector<SoilLayer>& profile, double freq_hz) {
  if (profile.empty()) throw std::invalid_argument("effectiveSoilTemperature: empty profile");

  // With alpha piecewise constant, the weight of layer i in Eq. (9) is
  // exp(-A_top) - exp(-A_bottom), A being the cumulative attenuation.
  const double k = 4.0 * kPi / wavelength(freq_hz);
  double tg = 0.0;
  double transmitted = 1.0;  // exp(-A) at the top of the current layer
  for (std::size_t i = 0; i < profile.size(); ++i) {
    const SoilLayer& layer = profile[i];
    if (i + 1 == profile.size()) {
      tg += layer.temp_k * transmitted;  // semi-infinite bottom layer absorbs the rest
      break;
    }
    const double alpha = k * layer.eps.imag() / (2.0 * std::sqrt(layer.eps.real()));  // Eq. (10)
    const double next = transmitted * std::exp(-alpha * layer.thickness_m);
    tg += layer.temp_k * (transmitted - next);
    transmitted = next;
  }
  return tg;
}

}  // namespace lmeb
