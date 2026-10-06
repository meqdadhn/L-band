// Unit tests for L-MEB (dependency-free).
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "lmeb/lmeb.hpp"

namespace {

int g_failures = 0;

#define CHECK(cond)                                                          \
  do {                                                                       \
    if (!(cond)) {                                                           \
      std::printf("  FAILED %s:%d: %s\n", __FILE__, __LINE__, #cond);        \
      ++g_failures;                                                          \
    }                                                                        \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                          \
  do {                                                                                 \
    const double va_ = (a), vb_ = (b);                                                 \
    if (!(std::abs(va_ - vb_) <= (tol))) {                                             \
      std::printf("  FAILED %s:%d: %s = %.6g, expected %.6g (tol %.3g)\n", __FILE__, \
                  __LINE__, #a, va_, vb_, static_cast<double>(tol));                   \
      ++g_failures;                                                                    \
    }                                                                                  \
  } while (0)

using namespace lmeb;

void testFreeWater() {
  // Pure water at 20 degC, 1.4 GHz: relaxation at ~17 GHz gives eps ~ 79.6 + 6.1i.
  const Complex e = freeWaterPermittivity(1.4e9, 20.0, 0.0);
  CHECK_NEAR(e.real(), 79.6, 0.5);
  CHECK_NEAR(e.imag(), 6.1, 0.3);
}

void testDobson() {
  const SoilTexture loam{0.3, 0.2, 1.3, 0.65};
  const Complex dry = dobsonPermittivity(0.02, 20.0, loam);
  const Complex wet = dobsonPermittivity(0.30, 20.0, loam);
  CHECK(dry.real() > 2.5 && dry.real() < 5.0);
  CHECK(wet.real() > 12.0 && wet.real() < 25.0);
  double prev = 0.0;
  for (double sm = 0.0; sm <= 0.5; sm += 0.05) {
    const double re = dobsonPermittivity(sm, 20.0, loam).real();
    CHECK(re > prev);
    prev = re;
  }
  // Frozen soil, Eq. (13).
  const Complex frozen = soilPermittivity(0.3, 260.0, loam);
  CHECK_NEAR(frozen.real(), 5.0, 1e-12);
  CHECK_NEAR(frozen.imag(), 0.5, 1e-12);
}

void testFresnel() {
  const Complex eps(20.0, 2.0);
  const DualPol r0 = fresnelReflectivity(eps, 0.0);
  CHECK_NEAR(r0.h, r0.v, 1e-12);  // no polarization difference at nadir
  // Real permittivity: r_V vanishes at the Brewster angle.
  const double brewster = std::atan(std::sqrt(9.0));
  CHECK_NEAR(fresnelReflectivity(Complex(9.0, 0.0), brewster).v, 0.0, 1e-12);
  // Nadir value |(1 - sqrt(eps)) / (1 + sqrt(eps))|^2.
  CHECK_NEAR(fresnelReflectivity(Complex(9.0, 0.0), 0.0).h, 0.25, 1e-12);
  // r_H increases with angle.
  CHECK(fresnelReflectivity(eps, deg2rad(50)).h > r0.h);
}

void testRoughness() {
  const DualPol fr{0.4, 0.2};
  // Eq. (8) with N_R = 0: uniform attenuation exp(-HR).
  const DualPol r = roughReflectivity(fr, {0.3, 0.0, 0.0, 0.0}, deg2rad(40));
  CHECK_NEAR(r.h, 0.4 * std::exp(-0.3), 1e-12);
  CHECK_NEAR(r.v, 0.2 * std::exp(-0.3), 1e-12);
  // N_R = 2 reproduces the classical cos^2 form.
  const double th = deg2rad(30);
  const DualPol r2 = roughReflectivity(fr, {0.3, 2.0, 2.0, 0.0}, th);
  CHECK_NEAR(r2.h, 0.4 * std::exp(-0.3 * std::pow(std::cos(th), 2)), 1e-12);
  // Q_R mixes polarizations.
  const DualPol rq = roughReflectivity(fr, {0.0, 0.0, 0.0, 0.5}, th);
  CHECK_NEAR(rq.h, 0.3, 1e-12);
  CHECK_NEAR(rq.v, 0.3, 1e-12);
  // Coherent component, Eq. (5): smooth surface unchanged.
  CHECK_NEAR(coherentReflectivity(fr, 0.0, th).h, 0.4, 1e-12);
}

void testEffectiveTemperature() {
  EffectiveTempParams p;  // w0 = 0.3, bw0 = 0.3
  CHECK_NEAR(effectiveTempCoefficient(0.3, p), 1.0, 1e-12);
  CHECK_NEAR(effectiveTempCoefficient(0.45, p), 1.0, 1e-12);  // clipped
  CHECK_NEAR(effectiveTempCoefficient(0.0, p), 0.0, 1e-12);
  CHECK_NEAR(effectiveSoilTemperature(0.3, 300.0, 290.0, p), 300.0, 1e-9);
  CHECK_NEAR(effectiveSoilTemperature(300.0, 290.0, 0.246), 292.46, 1e-9);

  // Eq. (9) with a uniform profile returns the profile temperature.
  const std::vector<SoilLayer> uniform = {
      {0.05, 290.0, {10.0, 2.0}}, {0.1, 290.0, {12.0, 2.5}}, {1.0, 290.0, {15.0, 3.0}}};
  CHECK_NEAR(effectiveSoilTemperature(uniform), 290.0, 1e-9);
  // A very lossy top layer dominates the emission.
  const std::vector<SoilLayer> lossy = {{0.5, 300.0, {30.0, 15.0}}, {1.0, 280.0, {5.0, 0.1}}};
  CHECK_NEAR(effectiveSoilTemperature(lossy), 300.0, 1e-3);
}

void testVegetation() {
  VegetationParams v;
  v.b_s1 = 0.05;
  v.b_s2 = -0.02;
  CHECK_NEAR(standingNadirOpticalDepth(2.0, v), 0.08, 1e-12);
  CHECK_NEAR(standingNadirOpticalDepth(0.1, v), 0.0, 1e-12);  // clipped at 0

  v.tt_h = 2.0;
  v.tt_v = 8.0;
  const DualPol t0 = standingOpticalDepth(0.2, 0.0, v);
  CHECK_NEAR(t0.h, 0.2, 1e-12);
  CHECK_NEAR(t0.v, 0.2, 1e-12);
  const DualPol t90 = standingOpticalDepth(0.2, kPi / 2, v);
  CHECK_NEAR(t90.h, 0.4, 1e-12);
  CHECK_NEAR(t90.v, 1.6, 1e-12);

  // Litter, Eqs. (25)-(27).
  v.a_l = 1.0;
  v.b_l = 0.2;
  v.litter_biomass = 0.5;
  CHECK_NEAR(litterMoistureContent(0.2, v), 0.4, 1e-12);
  CHECK_NEAR(litterMoistureContent(0.9, v), 0.8, 1e-12);  // clipped
  CHECK_NEAR(litterOpticalDepth(0.2, v), 0.24 * (0.4 / 0.6) * 0.5, 1e-12);

  // Composite temperature, Eqs. (16)-(18).
  CHECK_NEAR(compositeTemperature(300.0, 290.0, 0.0, 1.7), 290.0, 1e-12);
  CHECK_NEAR(compositeTemperature(300.0, 290.0, 5.0, 1.7), 300.0, 1e-12);  // A_t clipped
  const double at = 1.7 * (1.0 - std::exp(-0.3));
  CHECK_NEAR(compositeTemperature(300.0, 290.0, 0.3, 1.7), 290.0 + 10.0 * at, 1e-12);

  // Eq. (14) limits.
  CHECK_NEAR(tauOmegaBrightness(0.0, 1.0, 0.3, 280.0, 290.0), 0.7 * 290.0, 1e-9);  // bare
  CHECK_NEAR(tauOmegaBrightness(0.0, 0.0, 0.3, 280.0, 290.0), 280.0, 1e-9);        // opaque

  CHECK(interceptionFlag({280.0, 281.0}, VegetationParams{}));
  CHECK(!interceptionFlag({250.0, 280.0}, VegetationParams{}));
}

void testAtmosphere() {
  const DualPol tb = topOfAtmosphere({250.0, 270.0}, 0.01, 2.0, 0.0);
  CHECK_NEAR(tb.h, 250.0 * std::exp(-0.01) + 2.0, 1e-9);
  const DualPol agg = aggregatePixel({{0.25, {200.0, 220.0}}, {0.75, {240.0, 260.0}}});
  CHECK_NEAR(agg.h, 230.0, 1e-9);
  CHECK_NEAR(agg.v, 250.0, 1e-9);
}

void testForwardModel() {
  ModelConfig cfg = makeConfig(findCropPreset("corn"));
  SurfaceState st;
  st.sm = 0.25;
  st.tau_s_nad = 0.2;

  for (double th = 0.0; th <= 60.0; th += 10.0) {
    const ForwardResult r = simulate(cfg, st, deg2rad(th));
    CHECK(r.tb.h > 100.0 && r.tb.h < 320.0);
    CHECK(r.tb.v >= r.tb.h - 1e-9);  // V emission >= H for a land surface
  }
  // Drier soil is warmer, more vegetation is warmer.
  SurfaceState dry = st;
  dry.sm = 0.05;
  CHECK(simulate(cfg, dry, deg2rad(40)).tb.h > simulate(cfg, st, deg2rad(40)).tb.h);
  SurfaceState dense = st;
  dense.tau_s_nad = 0.6;
  CHECK(simulate(cfg, dense, deg2rad(40)).tb.h > simulate(cfg, st, deg2rad(40)).tb.h);

  // Bare soil: TB = (1 - r_G) T_G.
  ModelConfig bare_cfg = makeConfig(findCropPreset("bare"));
  SurfaceState bare = st;
  bare.tau_s_nad = 0.0;
  const ForwardResult b = simulate(bare_cfg, bare, deg2rad(30));
  CHECK_NEAR(b.tb.h, (1.0 - b.reflectivity.h) * b.t_ground, 1e-9);
  CHECK_NEAR(b.t_composite, b.t_ground, 1e-12);
}

void testRetrieval() {
  for (const char* crop : {"soybean", "wheat", "corn", "corn_rough"}) {
    const ModelConfig cfg = makeConfig(findCropPreset(crop));
    SurfaceState truth;
    truth.sm = 0.23;
    truth.tau_s_nad = 0.18;

    std::vector<Observation> obs;
    for (double th = 10.0; th <= 55.0; th += 5.0) {
      const DualPol tb = simulate(cfg, truth, deg2rad(th)).tb;
      obs.push_back({th, Pol::H, tb.h, 2.0});
      obs.push_back({th, Pol::V, tb.v, 2.0});
    }
    // Weak priors so that noise-free data drive the solution (SM and HR are
    // strongly correlated, so even sigma(SM) = 0.3 would bias rough soils).
    const std::vector<ParamPrior> priors = {
        {Param::SM, 0.05, 10.0},
        {Param::TauNad, 0.0, 10.0},
        {Param::HR, cfg.roughness.hr + 0.1, 10.0},
    };
    const RetrievalResult res = retrieve(cfg, truth, obs, priors);
    std::printf("  %-11s SM=%.4f tau=%.4f HR=%.4f rmse=%.2e it=%d\n", crop, res.value(Param::SM),
                res.value(Param::TauNad), res.value(Param::HR), res.rmse_tb, res.iterations);
    CHECK(res.converged);
    CHECK_NEAR(res.value(Param::SM), truth.sm, 2e-3);
    CHECK_NEAR(res.value(Param::TauNad), truth.tau_s_nad, 2e-3);
    CHECK_NEAR(res.value(Param::HR), cfg.roughness.hr, 1e-2);
    CHECK(res.rmse_tb < 0.05);
  }

  // '3-P' configuration (Table 3) on noisy synthetic data: HR constrained (sigma 0.1).
  {
    const ModelConfig cfg = makeConfig(findCropPreset("corn_rough"));
    SurfaceState truth;
    truth.sm = 0.27;
    truth.tau_s_nad = 0.25;
    std::vector<Observation> obs;
    const double noise[] = {1.1, -0.7, 0.4, -1.5, 0.9, 0.2, -0.3, 1.8, -1.0, 0.6};
    int k = 0;
    for (double th = 10.0; th <= 55.0; th += 5.0, ++k) {
      const DualPol tb = simulate(cfg, truth, deg2rad(th)).tb;
      obs.push_back({th, Pol::H, tb.h + noise[k % 10], 2.0});
      obs.push_back({th, Pol::V, tb.v - noise[(k + 3) % 10], 2.0});
    }
    const auto priors = threeParamPriors(cfg, truth, 0.0, cfg.roughness.hr);
    const RetrievalResult res = retrieve(cfg, truth, obs, priors);
    std::printf("  3-P noisy   SM=%.4f tau=%.4f HR=%.4f rmse=%.2f K\n", res.value(Param::SM),
                res.value(Param::TauNad), res.value(Param::HR), res.rmse_tb);
    CHECK(res.converged);
    CHECK_NEAR(res.value(Param::SM), truth.sm, 0.04);  // SMOS accuracy goal
    CHECK(res.rmse_tb < 2.0);
  }

  // The '3-P' priors of Table 3.
  const ModelConfig cfg = makeConfig(findCropPreset("soybean"));
  SurfaceState st;
  const auto priors = threeParamPriors(cfg, st, 0.0, 0.1);
  CHECK(priors.size() == 4);
  CHECK_NEAR(priors[0].initial, 0.05, 1e-12);
  CHECK_NEAR(priors[2].sigma, 0.1, 1e-12);
  CHECK_NEAR(priors[2].initial, effectiveSoilTemperature(st.sm, st.t_soil_surf, st.t_soil_depth,
                                                         cfg.eff_temp),
             1e-9);
}

}  // namespace

int main() {
  const std::vector<std::pair<std::string, std::function<void()>>> tests = {
      {"free water permittivity", testFreeWater},
      {"Dobson soil permittivity", testDobson},
      {"Fresnel reflectivity", testFresnel},
      {"roughness model", testRoughness},
      {"effective soil temperature", testEffectiveTemperature},
      {"vegetation tau-omega", testVegetation},
      {"atmosphere / aggregation", testAtmosphere},
      {"forward model", testForwardModel},
      {"3-P retrieval", testRetrieval},
  };
  for (const auto& [name, fn] : tests) {
    const int before = g_failures;
    std::printf("[ RUN  ] %s\n", name.c_str());
    fn();
    std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", name.c_str());
  }
  std::printf("%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
