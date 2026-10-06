#include "lmeb/presets.hpp"

#include <stdexcept>

namespace lmeb {

namespace {

VegetationParams veg(double b_s1, double b_s2, double tt_h, double tt_v, double omega) {
  VegetationParams v;
  v.b_s1 = b_s1;
  v.b_s2 = b_s2;
  v.tt_h = tt_h;
  v.tt_v = tt_v;
  v.omega_h = omega;
  v.omega_v = omega;
  return v;
}

// Soil textures of Table 2 (sand/clay as fractions, density in g/cm^3).
constexpr SoilTexture kAvignon{0.11, 0.27, 1.29, 0.65};    // PORTOS-91/-93
constexpr SoilTexture kEmirad{0.132, 0.328, 1.28, 0.65};   // EMIRAD-2001
constexpr SoilTexture kRebex{0.16, 0.29, 1.1, 0.65};       // REBEX
constexpr SoilTexture kBarc{0.67, 0.15, 1.2, 0.65};        // BARC

}  // namespace

const std::vector<CropPreset>& cropPresets() {
  // Roughness {HR, NRH, NRV, QR}; vegetation from Table 4 (tt, omega) and
  // Table 5 (b'_S, b''_S from the tau_NAD vs LAI regression).
  static const std::vector<CropPreset> presets = {
      {"bare", "Bare smooth soil (INRA Avignon)", {0.1, 0.0, 0.0, 0.0},
       veg(0.0, 0.0, 1.0, 1.0, 0.0), kAvignon},
      {"soybean", "Soybean, PORTOS-91", {0.1, 0.0, 0.0, 0.0},
       veg(0.090, -0.020, 1.0, 2.0, 0.0), kAvignon},
      {"wheat", "Wheat, PORTOS-93 (DoY < 161)", {0.1, 0.0, 0.0, 0.0},
       veg(0.034, -0.010, 1.0, 8.0, 0.0), kAvignon},
      {"corn", "Corn on smooth soil, EMIRAD-2001", {0.1, 0.0, 0.0, 0.0},
       veg(0.047, 0.0, 2.0, 1.0, 0.05), kEmirad},
      {"corn_rough", "Corn on rough soil, REBEX", {0.7, 0.5, -1.0, 0.0},
       veg(0.055, -0.018, 2.0, 1.0, 0.05), kRebex},
      {"soybean_barc", "Soybean, BARC (perpendicular)", {0.2, 0.0, -1.0, 0.0},
       veg(0.090, -0.020, 1.0, 1.0, 0.0), kBarc},
  };
  return presets;
}

const CropPreset& findCropPreset(const std::string& name) {
  for (const auto& p : cropPresets())
    if (p.name == name) return p;
  throw std::invalid_argument("unknown crop preset: " + name);
}

ModelConfig makeConfig(const CropPreset& preset) {
  ModelConfig cfg;
  cfg.texture = preset.texture;
  cfg.roughness = preset.roughness;
  cfg.vegetation = preset.vegetation;
  return cfg;
}

}  // namespace lmeb
