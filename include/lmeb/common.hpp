// Common types and physical constants for L-MEB.
#pragma once

#include <cmath>
#include <complex>

namespace lmeb {

using Complex = std::complex<double>;

constexpr double kPi = 3.14159265358979323846;
constexpr double kSpeedOfLight = 299792458.0;      // m/s
constexpr double kVacuumPermittivity = 8.854187817e-12;  // F/m
constexpr double kZeroCelsius = 273.15;            // K
constexpr double kSmosFrequencyHz = 1.4e9;         // L-band (1.4 GHz)

inline double deg2rad(double deg) { return deg * kPi / 180.0; }
inline double wavelength(double freq_hz) { return kSpeedOfLight / freq_hz; }

enum class Pol { H, V };

// A quantity defined at both horizontal and vertical polarization.
struct DualPol {
  double h = 0.0;
  double v = 0.0;

  double operator[](Pol p) const { return p == Pol::H ? h : v; }
  double& operator[](Pol p) { return p == Pol::H ? h : v; }
};

}  // namespace lmeb
