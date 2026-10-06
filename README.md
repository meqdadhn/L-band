# L-MEB — L-band Microwave Emission of the Biosphere (C++)

C++17 implementation of the L-MEB forward model and the '3-P' soil moisture
retrieval described in:

> Wigneron, J.-P. et al. (2007). *L-band Microwave Emission of the Biosphere
> (L-MEB) Model: Description and calibration against experimental data sets
> over crop fields.* Remote Sensing of Environment 107, 639–655.
> doi:10.1016/j.rse.2006.10.014

No external dependencies.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Options: `-DLMEB_BUILD_APPS=OFF`, `-DLMEB_BUILD_TESTS=OFF`.

## Tools

```bash
# TB vs incidence angle for a crop preset
./build/lmeb_forward --crop corn --sm 0.25 --lai 3 --tsurf 296 --tdepth 290 --tc 294

# '3-P' retrieval (SM, tau_NAD, HR) on synthetic noisy data
./build/lmeb_retrieve --crop wheat --sm 0.18 --tau 0.12 --noise 2

# ... or on measured data: CSV lines "theta_deg,pol,tb[,sigma]"
./build/lmeb_retrieve --crop soybean --obs obs.csv --tsurf 296 --tdepth 291 --tc 295
```

Run either tool with `--help` for all options.

## Running the τ–ω model

### What it computes

L-MEB treats the vegetation as one layer over the soil and computes the
surface brightness temperature at each polarization P (H or V) with the
τ–ω model (Eq. 14):

```
TB_P = (1−ω_P)(1−γ_P)(1 + γ_P r_GP) T_C  +  (1 − r_GP) γ_P T_G
        └── canopy emission, direct + ──┘   └─ soil emission ─┘
            reflected by the soil             through the canopy
```

| Symbol | Meaning | Where it comes from |
|---|---|---|
| τ_P | vegetation optical depth | τ_NAD (from `--tau`, or from LAI via Eq. 22), made angle- and polarization-dependent by tt_H, tt_V (Eqs. 23–24) |
| γ_P | canopy transmissivity, exp(−τ_P / cos θ) | Eq. 15 |
| ω_P | single-scattering albedo | crop preset or `--omega` |
| r_GP | rough-soil reflectivity | soil moisture → permittivity → Fresnel → roughness H_R, N_RP |
| T_G | effective soil temperature | `--tsurf`, `--tdepth` and soil moisture (Eqs. 11–12) |
| T_C | canopy temperature | `--tc` |

L-MEB replaces both T_C and T_G with the composite temperature
T_GC = A_t T_C + (1 − A_t) T_G (Eqs. 16–17).

Two limits help as a sanity check:
- **τ = 0 (bare soil):** γ = 1, so TB = (1 − r_G) T_G.
- **τ → ∞ (opaque canopy):** γ = 0, so TB = (1 − ω) T_C.

### Where it is in the code

| Step | Function | File |
|---|---|---|
| Eq. 14, the τ–ω equation | `tauOmegaBrightness` | `src/vegetation.cpp` |
| Eq. 15, γ_P | `vegetationAttenuation` | `src/vegetation.cpp` |
| Eqs. 23–24, τ_P(θ) | `standingOpticalDepth` | `src/vegetation.cpp` |
| Eq. 22, τ_NAD from LAI | `standingNadirOpticalDepth` | `src/vegetation.cpp` |
| Eqs. 16–17, T_GC | `compositeTemperature` | `src/vegetation.cpp` |
| Soil + vegetation combined | `simulate` | `src/model.cpp` |

### Step 1: build

```bash
cd ~/Dev/my_repos/lmeb
cmake -S . -B build
cmake --build build -j
```

### Step 2: run it for one canopy over a range of angles

```bash
./build/lmeb_forward --crop soybean --sm 0.3 --tau 0.3 --theta-step 10
```

```
# L-MEB forward model: soybean (Soybean, PORTOS-91)
# SM=0.300  eps_G=14.49+2.95i  T_G=295.00 K  T_GC=294.12 K  tau_NAD=0.300
# HR=0.10 NRH=0.00 NRV=0.00 tt_H=1.00 tt_V=2.00 omega=0.000/0.000
   theta       TBH       TBV      rGH      rGV     tauH     tauV       PR
     0.0   243.433   243.433   0.3140   0.3140   0.3000   0.3000   0.0000
    10.0   243.104   245.596   0.3190   0.3090   0.3000   0.3090   0.0051
    20.0   242.200   251.784   0.3343   0.2937   0.3000   0.3351   0.0194
    30.0   241.024   261.100   0.3609   0.2669   0.3000   0.3750   0.0400
    40.0   240.265   272.074   0.4007   0.2267   0.3000   0.4240   0.0621
    50.0   241.356   282.692   0.4562   0.1709   0.3000   0.4760   0.0789
    60.0   247.109   290.560   0.5307   0.0988   0.3000   0.5250   0.0808
```

The header shows the intermediate values: the soil permittivity ε_G, the
effective temperatures T_G and T_GC, and τ_NAD. Each row gives, for one
incidence angle:

| Column | Meaning |
|---|---|
| `TBH`, `TBV` | brightness temperature in K, the output of the τ–ω equation |
| `rGH`, `rGV` | rough-soil reflectivity r_GP |
| `tauH`, `tauV` | optical depth τ_P at that angle. Here only τ_V grows, because soybean has tt_V = 2 and tt_H = 1. |
| `PR` | polarization ratio (TBV − TBH) / (TBV + TBH) |

### Step 3: vary the vegetation

Sweeping τ_NAD shows the τ–ω model moving from soil-dominated to
canopy-dominated emission (soybean, SM = 0.3, θ = 40°):

```bash
for t in 0 0.2 0.5 1 3; do
  ./build/lmeb_forward --crop soybean --sm 0.3 --tau $t --theta-max 40 --theta-step 40 | tail -1
done
```

| τ_NAD | TBH (K) | TBV (K) | PR |
|---|---|---|---|
| 0 (bare) | 176.8 | 228.1 | 0.127 |
| 0.2 | 224.4 | 262.5 | 0.078 |
| 0.5 | 261.8 | 283.1 | 0.039 |
| 1 | 284.4 | 291.3 | 0.012 |
| 3 | 293.0 | 293.0 | 0.000 |

As the canopy thickens, TB rises toward the canopy temperature, the H/V
difference disappears, and sensitivity to soil moisture is lost.

### Inputs you can change

| Option | Parameter | Default |
|---|---|---|
| `--crop` | preset for H_R, N_RP, tt_P, ω, b′_S, b″_S and soil texture | `corn` |
| `--sm` | soil moisture (m³/m³) | 0.2 |
| `--tau` | τ_NAD, nadir optical depth | (from LAI) |
| `--lai` | LAI; τ_NAD = b′_S·LAI + b″_S, used when `--tau` is not given | 2 |
| `--omega` | ω (both polarizations) | preset |
| `--tth`, `--ttv` | angular dependence of τ at H and V | preset |
| `--hr`, `--nrh`, `--nrv` | soil roughness | preset |
| `--tsurf`, `--tdepth`, `--tc` | soil surface, deep soil and canopy temperature (K) | 295, 290, 293 |
| `--sand`, `--clay`, `--density` | soil texture (fractions, g/cm³) | preset |
| `--tau-atm`, `--tb-sky` | apply the atmosphere correction (Eq. 1) | off |
| `--theta-max`, `--theta-step` | angle sweep (degrees) | 60, 5 |

Examples:

```bash
# Corn canopy from LAI, with scattering (ω = 0.05) and tt_H = 2
./build/lmeb_forward --crop corn --sm 0.25 --lai 3

# Wheat: strong increase of τ_V with angle (tt_V = 8)
./build/lmeb_forward --crop wheat --sm 0.2 --tau 0.2

# Isotropic canopy with stronger scattering
./build/lmeb_forward --crop soybean --sm 0.2 --tau 0.4 --tth 1 --ttv 1 --omega 0.1
```

### Calling the τ–ω equation directly from C++

To bypass the soil model and run the τ–ω equation with your own reflectivity
and temperatures:

```cpp
#include <lmeb/lmeb.hpp>
#include <cstdio>

int main() {
  const double theta = lmeb::deg2rad(40.0);

  lmeb::VegetationParams veg;
  veg.tt_h = 1.0;
  veg.tt_v = 2.0;

  // Optical depth at this angle, then canopy transmissivity (Eqs. 23-24, 15)
  lmeb::DualPol tau = lmeb::standingOpticalDepth(/*tau_nad=*/0.3, theta, veg);
  lmeb::DualPol gamma = lmeb::vegetationAttenuation(tau, theta);

  // Eq. 14 with given soil reflectivity and temperatures
  double tbh = lmeb::tauOmegaBrightness(/*omega=*/0.0, gamma.h, /*r_soil=*/0.40,
                                        /*T_C=*/293.0, /*T_G=*/295.0);
  double tbv = lmeb::tauOmegaBrightness(0.0, gamma.v, 0.23, 293.0, 295.0);
  std::printf("TBH = %.2f K, TBV = %.2f K\n", tbh, tbv);
}
```

## Library usage

```cpp
#include <lmeb/lmeb.hpp>

lmeb::ModelConfig cfg = lmeb::makeConfig(lmeb::findCropPreset("soybean"));
lmeb::SurfaceState st;
st.sm = 0.22;
st.tau_s_nad = lmeb::standingNadirOpticalDepth(/*LAI=*/2.5, cfg.vegetation);
lmeb::DualPol tb = lmeb::simulate(cfg, st, lmeb::deg2rad(40.0)).tb;
```

Link with `lmeb::lmeb` (`add_subdirectory(lmeb)` or `find_package(lmeb)` after install).

## Equation map

| Paper | Code |
|---|---|
| Eq. 1 top-of-atmosphere TB | `topOfAtmosphere` (atmosphere.hpp) |
| Mixed-pixel linear combination | `aggregatePixel` |
| Eq. 2–3 emissivity, Fresnel | `fresnelReflectivity` (soil.hpp) |
| Eq. 5 coherent reflectivity | `coherentReflectivity` |
| Eq. 6–8 roughness (H_R, N_RP, Q_R) | `roughReflectivity` |
| Eq. 9–10 RT effective temperature | `effectiveSoilTemperature(profile)` |
| Eq. 11–12 T_G, C_t = (SM/w0)^bw0 | `effectiveSoilTemperature(sm, ...)` |
| Dobson et al. (1985) permittivity | `dobsonPermittivity` (dielectric.hpp) |
| Eq. 13 frozen soil ε = 5 + 0.5i | `soilPermittivity` (T < 273.15 K) |
| Eq. 14 τ–ω model | `tauOmegaBrightness` (vegetation.hpp) |
| Eq. 15 γ_P | `vegetationAttenuation` |
| Eq. 16–18 composite T_GC | `compositeTemperature` |
| Eq. 19 τ_P = τ_SP + τ_L + τ_IP | `simulate` (model.cpp) |
| Eq. 22 τ_S_NAD = b'_S LAI + b''_S | `standingNadirOpticalDepth` |
| Eq. 23–24 tt_H, tt_V | `standingOpticalDepth` |
| Eq. 25–27 litter | `litterOpticalDepth`, `litterMoistureContent` |
| PR interception flag (PR_I = 0.02) | `interceptionFlag` |
| Eq. 28 cost function, Table 3 '3-P' | `retrieve`, `threeParamPriors` (retrieval.hpp) |
| Tables 2, 4, 5 calibrations | `cropPresets` (presets.cpp) |

The retrieval minimizes Eq. 28 with Levenberg–Marquardt (numerical Jacobian,
box bounds on each parameter) and reports posterior standard deviations from
(JᵀJ)⁻¹. Any subset of SM, τ_NAD, T_GC, H_R, N_RH, N_RV, ω_H, ω_V, tt_H, tt_V
can be retrieved; the others stay fixed at their `ModelConfig`/`SurfaceState` values.

## Crop presets

| Name | Source | H_R | N_RH | N_RV | tt_H | tt_V | ω | b'_S | b''_S |
|---|---|---|---|---|---|---|---|---|---|
| `bare` | Avignon soil | 0.1 | 0 | 0 | 1 | 1 | 0 | – | – |
| `soybean` | PORTOS-91 | 0.1 | 0 | 0 | 1 | 2 | 0 | 0.090 | −0.020 |
| `wheat` | PORTOS-93 | 0.1 | 0 | 0 | 1 | 8 | 0 | 0.034 | −0.010 |
| `corn` | EMIRAD-2001 | 0.1 | 0 | 0 | 2 | 1 | 0.05 | 0.047 | 0.0 |
| `corn_rough` | REBEX | 0.7 | 0.5 | −1 | 2 | 1 | 0.05 | 0.055 | −0.018 |
| `soybean_barc` | BARC ⊥ | 0.2 | 0 | −1 | 1 | 1 | 0 | 0.090* | −0.020* |

\* LAI was not available for BARC; the PORTOS-91 soybean LAI regression is reused.
Soil textures come from Table 2 (salinity 0.65 ‰). Defaults: w0 = bw0 = 0.3, B_t = 1.7, c_L = 0.24.

## Choices not specified in the paper

- **Dobson model details:** Dobson et al. (1985) with α = 0.65, ρ_s = 2.66 g/cm³ and
  σ_eff = −1.645 + 1.939ρ_b − 2.256S + 1.594C (1.4–18 GHz form). Free water uses the
  Debye model with Ulaby's temperature fits and Klein–Swift salinity corrections.
- **Soil permittivity temperature:** computed at the surface soil temperature;
  soil is treated as frozen when that temperature is below 0 °C.
- **Composite temperature:** `simulate` uses T_GC for both T_C and T_G in Eq. 14
  (set `use_composite_temperature = false` to use them separately). τ_NAD in Eq. 17 is
  the total nadir optical depth (standing vegetation + litter + interception).
- **Not implemented:** Mätzler's (1998) dry-sand permittivity, Pellarin et al.'s
  (2003a) parameterization of τ_ATM and TB_SKY_U (these are inputs here), the antenna
  pattern weighting, the forest and other generic classes, and the topography index.

The SM–H_R trade-off is close to degenerate on rough soils. That is why '3-P'
constrains H_R (σ = 0.1). With looser priors, expect the retrieval to move along that ridge.
