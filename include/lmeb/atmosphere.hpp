// Atmospheric correction and pixel aggregation (Section 2.1).
#pragma once

#include <utility>
#include <vector>

#include "lmeb/common.hpp"

namespace lmeb {

// TB_P_TOA = TB_P exp(-tau_ATM / cos(theta)) + TB_SKY_U, Eq. (1).
DualPol topOfAtmosphere(const DualPol& tb_surface, double tau_atm, double tb_sky_up,
                        double theta_rad);

// Brightness temperature of a mixed pixel: linear combination of the TB of each
// cover type weighted by its cover fraction (fractions are normalized to sum 1).
DualPol aggregatePixel(const std::vector<std::pair<double, DualPol>>& fraction_and_tb);

}  // namespace lmeb
