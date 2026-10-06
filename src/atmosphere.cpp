#include "lmeb/atmosphere.hpp"

#include <stdexcept>

namespace lmeb {

DualPol topOfAtmosphere(const DualPol& tb, double tau_atm, double tb_sky_up, double theta_rad) {
  const double t = std::exp(-tau_atm / std::cos(theta_rad));
  return {tb.h * t + tb_sky_up, tb.v * t + tb_sky_up};
}

DualPol aggregatePixel(const std::vector<std::pair<double, DualPol>>& fraction_and_tb) {
  double total = 0.0;
  DualPol out;
  for (const auto& [frac, tb] : fraction_and_tb) {
    total += frac;
    out.h += frac * tb.h;
    out.v += frac * tb.v;
  }
  if (total <= 0.0) throw std::invalid_argument("aggregatePixel: cover fractions sum to zero");
  out.h /= total;
  out.v /= total;
  return out;
}

}  // namespace lmeb
