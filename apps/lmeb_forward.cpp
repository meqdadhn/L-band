// Forward simulation of L-band brightness temperatures with L-MEB.
//
// Example:
//   lmeb_forward --crop corn --sm 0.25 --lai 3 --tsurf 296 --tdepth 290 --tc 294
#include <cstdio>
#include <iostream>

#include "cli.hpp"
#include "lmeb/lmeb.hpp"

namespace {

void usage() {
  std::cout <<
      "Usage: lmeb_forward [options]\n"
      "  --crop NAME       bare|soybean|wheat|corn|corn_rough|soybean_barc (default corn)\n"
      "  --sm X            soil moisture, m3/m3 (default 0.2)\n"
      "  --lai X           leaf area index; tau_S_NAD = b'_S LAI + b''_S (default 2)\n"
      "  --tau X           nadir optical depth (overrides --lai)\n"
      "  --tsurf K         surface soil temperature (default 295)\n"
      "  --tdepth K        deep soil temperature (default 290)\n"
      "  --tc K            canopy temperature (default 293)\n"
      "  --hr X --nrh X --nrv X --omega X --tth X --ttv X   override preset parameters\n"
      "  --sand F --clay F --density G                     override soil texture\n"
      "  --tau-atm X --tb-sky K   top-of-atmosphere TB, Eq. (1) (default: off)\n"
      "  --theta-max DEG --theta-step DEG                  angle sweep (default 0..60 step 5)\n";
}

}  // namespace

int main(int argc, char** argv) {
  cli::Args args(argc, argv);
  if (args.has("help")) {
    usage();
    return 0;
  }

  lmeb::CropPreset preset;
  try {
    preset = lmeb::findCropPreset(args.str("crop", "corn"));
  } catch (const std::exception& e) {
    cli::Args::fail(e.what());
  }
  lmeb::ModelConfig cfg = lmeb::makeConfig(preset);
  cfg.roughness.hr = args.num("hr", cfg.roughness.hr);
  cfg.roughness.nr_h = args.num("nrh", cfg.roughness.nr_h);
  cfg.roughness.nr_v = args.num("nrv", cfg.roughness.nr_v);
  cfg.vegetation.tt_h = args.num("tth", cfg.vegetation.tt_h);
  cfg.vegetation.tt_v = args.num("ttv", cfg.vegetation.tt_v);
  cfg.vegetation.omega_h = cfg.vegetation.omega_v =
      args.num("omega", cfg.vegetation.omega_h);
  cfg.texture.sand = args.num("sand", cfg.texture.sand);
  cfg.texture.clay = args.num("clay", cfg.texture.clay);
  cfg.texture.bulk_density = args.num("density", cfg.texture.bulk_density);

  lmeb::SurfaceState st;
  st.sm = args.num("sm", 0.2);
  st.t_soil_surf = args.num("tsurf", 295.0);
  st.t_soil_depth = args.num("tdepth", 290.0);
  st.t_canopy = args.num("tc", 293.0);
  st.tau_s_nad = args.has("tau")
                     ? args.num("tau", 0.0)
                     : lmeb::standingNadirOpticalDepth(args.num("lai", 2.0), cfg.vegetation);

  const bool toa = args.has("tau-atm");
  const double tau_atm = args.num("tau-atm", 0.0);
  const double tb_sky = args.num("tb-sky", 0.0);
  const double theta_max = args.num("theta-max", 60.0);
  const double theta_step = args.num("theta-step", 5.0);
  if (theta_step <= 0.0) cli::Args::fail("--theta-step must be > 0");

  const lmeb::ForwardResult nadir = lmeb::simulate(cfg, st, 0.0);
  std::printf("# L-MEB forward model: %s (%s)\n", preset.name.c_str(), preset.description.c_str());
  std::printf("# SM=%.3f  eps_G=%.2f+%.2fi  T_G=%.2f K  T_GC=%.2f K  tau_NAD=%.3f\n", st.sm,
              nadir.eps_soil.real(), nadir.eps_soil.imag(), nadir.t_ground, nadir.t_composite,
              nadir.tau_nad);
  std::printf("# HR=%.2f NRH=%.2f NRV=%.2f tt_H=%.2f tt_V=%.2f omega=%.3f/%.3f\n",
              cfg.roughness.hr, cfg.roughness.nr_h, cfg.roughness.nr_v, cfg.vegetation.tt_h,
              cfg.vegetation.tt_v, cfg.vegetation.omega_h, cfg.vegetation.omega_v);
  std::printf("%8s %9s %9s %8s %8s %8s %8s %8s\n", "theta", "TBH", "TBV", "rGH", "rGV", "tauH",
              "tauV", "PR");
  for (double th = 0.0; th <= theta_max + 1e-9; th += theta_step) {
    const lmeb::ForwardResult r = lmeb::simulate(cfg, st, lmeb::deg2rad(th));
    const lmeb::DualPol tb =
        toa ? lmeb::topOfAtmosphere(r.tb, tau_atm, tb_sky, lmeb::deg2rad(th)) : r.tb;
    std::printf("%8.1f %9.3f %9.3f %8.4f %8.4f %8.4f %8.4f %8.4f\n", th, tb.h, tb.v,
                r.reflectivity.h, r.reflectivity.v, r.tau.h, r.tau.v,
                lmeb::polarizationRatio(tb));
  }
  return 0;
}
