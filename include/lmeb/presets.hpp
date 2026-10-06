// Calibrated parameters for crops (Tables 2, 4 and 5 of Wigneron et al. 2007).
#pragma once

#include <string>
#include <vector>

#include "lmeb/model.hpp"

namespace lmeb {

struct CropPreset {
  std::string name;
  std::string description;
  RoughnessParams roughness;
  VegetationParams vegetation;
  SoilTexture texture;
};

// Available presets: bare, soybean, wheat, corn, corn_rough, soybean_barc.
const std::vector<CropPreset>& cropPresets();

// Throws std::invalid_argument for unknown names.
const CropPreset& findCropPreset(const std::string& name);

// ModelConfig built from a preset (L-band, composite temperature).
ModelConfig makeConfig(const CropPreset& preset);

}  // namespace lmeb
