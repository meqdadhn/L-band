// Soil permittivity (Section 2.2.3).
//
// Unfrozen soil: semi-empirical mixing model of Dobson et al. (1985).
// Frozen soil:   constant eps_GF = 5 + 0.5i (Eq. 13).
//
// Permittivities use the convention eps = eps' + i*eps'' (eps'' >= 0).
#pragma once

#include "lmeb/common.hpp"

namespace lmeb {

struct SoilTexture {
  double sand = 0.30;          // sand mass fraction [0-1]
  double clay = 0.20;          // clay mass fraction [0-1]
  double bulk_density = 1.30;  // g/cm^3
  double salinity = 0.65;      // per mil (paper uses 0.65 for all data sets)
};

// Permittivity of eps_GF of frozen soil, Eq. (13).
constexpr Complex kFrozenSoilPermittivity{5.0, 0.5};

// Permittivity of (slightly saline) free water, Debye relaxation with the
// Stogryn/Klein-Swift temperature and salinity dependencies.
Complex freeWaterPermittivity(double freq_hz, double temp_c, double salinity_ppt);

// Dobson et al. (1985) mixing model. sm is volumetric soil moisture (m3/m3).
Complex dobsonPermittivity(double sm, double temp_c, const SoilTexture& texture,
                           double freq_hz = kSmosFrequencyHz);

// Soil permittivity eps_G: Dobson model above 0 degC, Eq. (13) below.
Complex soilPermittivity(double sm, double temp_k, const SoilTexture& texture,
                         double freq_hz = kSmosFrequencyHz);

}  // namespace lmeb
