# Plume-surface interaction validation

This study evaluates SpaceAGORA's plume-surface interaction model
(`src/dynamics/coupled/force_torque_models/plume_surface_interaction.jl`,
documented in [`docs/src/user/lunar_landing.md`](../../../docs/src/user/lunar_landing.md))
against every published plume-surface measurement that could be sourced
precisely, and records the gaps.

It is modeled on [`benchmarks/studies/telemetry_validation/`](../telemetry_validation/README.md):
one frozen manifest per reference case under `manifests/`, a runner that writes
a comparison table to `results/` (gitignored), and `common.jl` for the shared
parts. Unlike that study it runs no simulation — every case is a direct
evaluation of the model's public functions (`plume_surface_footprint`,
`plume_quantities`, `plume_erosion_onset_height`, `plume_ground_effect_force`,
`plume_ejecta_summary`) at the case's conditions, with two cases integrating a
published descent profile in time (the total eroded mass and the crater depth,
through `model_profile_mass` and `model_profile_crater` in `common.jl`). The
whole thing takes a couple of seconds after Julia loads, needs no SPICE, no GRAM
and no results file.

## Commands

```bash
# Every case on the default analytic plume field (nothing gates)
julia --project=. benchmarks/studies/psi_validation/run_validation.jl

# The same model on the tabulated source-flow plume of the LM descent engine
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --field=table

# Only the cases that may gate, and fail the run if one is outside tolerance
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --cases=gates --enforce=true

# One case
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --cases=pfgt1_erosion_threshold_shear
```

`--field` takes `analytic` (the default `PlumeAnalyticField`), `table`
(`data/psi/apollo_lmde.json`) or a path to another plume field table. Everything
else about the model is the shipped default: the study reports the model as it
is, never as it could be tuned.

Outputs (gitignored):

```
results/<hostname>/<analytic|table>/
  psi_validation_comparison.csv   # one row per reference value: model, reference, ratio, verdict
  psi_validation_cases.csv        # one row per case: status, gating, tolerance, full citation
```

The run ends with a findings block computed from the rows it just produced: the
model/reference span of every case, whether the model's trend with height has
the same sign as the reference's, where each series peaks, and what is on record
as unmeasurable. The prose below is written from that block and the table above
it, not from expectation.

Four of the cases also run as PR-tier regression tests in
`test/unit/dynamics/psi_reference_tests.jl`, which loads these same manifests
rather than repeating their numbers, so a reference value and its source can
never drift apart between the study and the test.

## Rules this study follows

1. **No number without a source.** Every reference value names its authors,
   title, venue, year and the table, figure, equation or section it sits on, in
   the manifest's `[source]` block. Values the harness derives from a published
   table (for example a product of two published columns) are flagged
   `derived = true` with the arithmetic written out.
2. **A case that cannot be sourced does not gate.** Each manifest carries a
   `case.status`, and only `sourced` cases may set `gate_eligible = true`; the
   loader refuses any other combination. The five non-gating statuses are
   `cross_regime` (a real measurement from a different body or flow regime),
   `calibration_target` (the observable a model parameter was fitted to, so
   agreement is circular), `different_quantity` (the model produces something
   close but not a definition of the same thing, so no tolerance can be argued),
   `not_modeled` (the model produces nothing comparable) and `unsourced` (no
   published value was found at all).
3. **The tolerance comes from the reference, not from the model.** Each
   manifest's `tolerance.reason` argues the band from the reference's own
   scatter, stated uncertainty, or the spread between independent published
   estimates of the same quantity. Nothing here was widened to make a case pass.
4. **The model was not touched.** These numbers are today's model with today's
   defaults. Two cases fail; the failures are reported, not fixed.

## Reference cases

| Case | Quantity | Source | Gate-eligible |
|---|---|---|---|
| `apollo12_erosion_rate` | mass erosion rate vs. engine height, 11 altitudes | Lane & Metzger 2015, Table 2 + Eq. (11) | **yes** |
| `apollo12_scour_radius` | radius of the eroding region vs. height, 11 altitudes | Lane & Metzger 2015, Table 2 column `a0` | **yes** |
| `apollo_total_eroded_mass` | total soil moved in one landing | Lane & Metzger 2015, Table 3 (+ 4 other published estimates) | **yes** |
| `apollo_lm_surface_shear` | wall shear stress vs. height | Lane & Metzger 2015, Eq. (14) / Fig. 11 | **yes** |
| `pfgt1_erosion_threshold_shear` | shear at which erosion starts | Stubbs & Mehta 2026 (NTRS 20250011216), Fig. 14 | **yes** |
| `surveyor3_ejecta_speed` | grain speed | Immer et al. 2008, §1.1 (4 primary estimates) | **yes** |
| `apollo_scour_depth` | depth of the toroidal scour | Metzger 2024 part 2, §5 | **yes** (new) |
| `apollo16_erosion_onset_height` | height at which dust starts blowing | Metzger 2024 part 2, §2 | no — different quantity |
| `viking_deep_cratering_pressure` | surface pressure at which deep cratering starts | Metzger 2016, "Deep cratering on Mars" | no — cross-regime |
| `apollo_dust_ejection_angle` | ejection angle above horizontal | Immer et al. 2008, Table 1 | no — calibration target |
| `ground_effect_augmentation` | thrust augmentation near the ground | none found | no — **unsourced** |

Full citations, conditions, uncertainties and tolerance arguments are in the
manifests; the run prints the citation line for every case it evaluates.

### The sources

- **J. E. Lane and P. T. Metzger**, "Estimation of Apollo lunar dust transport
  using optical extinction measurements", *Acta Geophysica* **63**(2), 568–599
  (2015), doi:10.1515/acgeo-2015-0005; open preprint
  [arXiv:1503.00154](https://arxiv.org/abs/1503.00154). The single richest
  source here: an inversion of the optical extinction of the blowing dust in the
  Apollo 12 descent film into an erosion rate at eleven altitudes (Table 2), a
  total eroded mass (Table 3), and a CFD surface-shear fit
  `τ(h) = 6.21 exp(−0.123 h)` Pa (Eq. 14, Fig. 11). Its `h` is defined in the
  nomenclature as the *nozzle opening distance from the surface*, which is the
  height the model's `plume_quantities` takes.
- **P. T. Metzger**, "Erosion rate of lunar soil under a landing rocket, part 2:
  Benchmarking and predictions", *Icarus* **417**, 116135 (2024),
  doi:10.1016/j.icarus.2024.116135; open preprint
  [arXiv:2403.18584](https://arxiv.org/abs/2403.18584). Gives the Apollo 16
  onset of visible dust at 31.5 m, the 3–12 cm post-landing scour depth (which
  it attributes to Metzger, Smith & Lane 2011, *JGR* **116**, E06005), and a
  modern total-erosion estimate of 11–26 t per Apollo landing, four to ten times
  the earlier numbers.
- **D. C. Stubbs and M. Mehta**, "A Data-Derived Scaling Approach for
  Plume-Surface Interaction Crater Formation", AIAA SciTech Forum 2026, NASA
  Technical Reports Server [20250011216](https://ntrs.nasa.gov/citations/20250011216)
  (open full text). The only sourced measurement of an erosion *threshold*
  shear stress: Loci/CHEM CFD shear evaluated at the optically measured crater
  edge in NASA MSFC's Physics Focused Ground Test 1 vacuum-chamber campaign,
  averaging 0.25 Pa on mono-disperse silica sand and 0.47 Pa on sieved BP-1.
- **C. Immer, J. Lane, P. Metzger and S. Clements**, "Apollo Video
  Photogrammetry Estimation of Plume Impingement Effects", Earth & Space 2008
  (11th ASCE Aerospace Division International Conference, Long Beach); open
  preprint [arXiv:2104.07669](https://arxiv.org/abs/2104.07669). Table 1 gives
  the per-mission dust ejection angle (Apollo 11 2.6°, 14 2.4°, 15 8.1°, 16
  1.4°, 17 2.0°); §1.1 compiles the four Surveyor III ejecta-speed estimates.
- **P. T. Metzger**, "ISRU Implications for Lunar and Martian Plume Effects",
  AIAA SciTech 2016; open preprint
  [arXiv:2104.06248](https://arxiv.org/abs/2104.06248). The 0.3 kPa Viking
  empirical deep-cratering limit, and the statement that the near-ground base
  pressure rise "has not been characterized" — the evidence for the
  ground-effect case being unsourced.
- **P. T. Metzger, J. E. Lane, C. D. Immer and X. Li**, "Scaling of Erosion Rate
  in Subsonic Jet Experiments and Apollo Lunar Module Landings", Earth & Space
  2010; open preprint [arXiv:2104.10038](https://arxiv.org/abs/2104.10038).
  Context for the shear case: its RANS CFD of the LM plume integrates to
  17.4 mN at 20 m and 35.9 mN at 10 m engine height, and it states plainly that
  the smooth-wall CFD understates the real shear, that Roberts' rough-wall
  approximation (shear ≈ dynamic pressure) is orders of magnitude higher, and
  that the result "is not to be considered a physical erosion law".
- **R. W. Orloff**, *Apollo by the Numbers: A Statistical Reference*, NASA
  SP-2000-4029, "Selected Mission Weights (lbs)". The landed LM masses the
  near-hover thrusts are derived from: Apollo 11 16,153.2 lb, Apollo 12
  16,564.2 lb, Apollo 16 18,208 lb.
- **E. M. Jones and K. Glover** (eds.), *Apollo 11 Lunar Surface Journal*, NASA,
  landing page. Supporting context for the onset case: Armstrong's technical
  debrief, "I first noticed that we were, in fact, disturbing the dust on the
  surface when we were something less than 100 feet", and Aldrin at 102:45:17,
  "40 feet, down 2 1/2. Picking up some dust."

### Cases that were looked for and are not here

- **NASA Langley 60-foot vacuum sphere (HLS PSI Ground Test).** The campaign is
  running and will instrument a flat plate with pressure, shear and heat-flux
  sensors, which would make it the best possible reference for this model. What
  is published openly so far is flow visualization — laser Rayleigh scattering
  giving the stagnation shock's standoff height and curvature — and the model
  produces no such quantity. The initial-measurements paper (AIAA SciTech 2026)
  is behind the AIAA paywall. Nothing from it is used, and no number from it is
  quoted anywhere in this study.
- **Raw MSFC PFGT-1 crater depths and widths.** The open campaign overview
  (Eberhart, West & Korzun, NTRS 20210025519) is qualitative; the crater depth
  and width time histories appear only inside figures of follow-on papers. The
  one quantitative result that could be read from the open literature is the
  crater-edge threshold shear, which is the `pfgt1_erosion_threshold_shear`
  case. (The model does now produce a crater depth; it is scored against
  Metzger's Apollo interpretation in `apollo_scour_depth`, the only sourced
  crater depth this study could find.)
- **LROC blast-zone and regolith-brightening radii at the Apollo and Chang'e-3
  sites.** The primary papers (Clegg et al., *Icarus* 2014; Clegg-Watkins et
  al., *Icarus* 2016) are paywalled, the secondary summaries that could be read
  disagreed with each other by nearly an order of magnitude on the blast-zone
  area, and the brightening halo is a far-field deposition signature rather than
  the locally eroded annulus the model computes. Excluded as unsourced rather
  than guessed.
- **Phoenix on Mars** (Mehta et al., *Icarus* **211**, 172–194, 2011). Pulsed
  monopropellant thrusters in a Martian atmosphere producing explosive erosion —
  a regime outside the model's steady single-jet closure. Not used.

## Baseline: the single fitted law (commit `5c4ed090`)

Recorded on `space-falcon-1` from the model as it stood before the erosion
regimes were wired in: one Roberts shear-excess law with two calibrated
constants (`threshold_shear_pa` 0.15 fitted to put the onset at 31 m,
`erosion_efficiency` 10 standing in for the saltation cascade), on the Gaussian
`PlumeAnalyticField`. **Four of six gate-eligible cases passed; two failed.**
This is the reference point everything below has to beat.

| Case | Model | Reference | Model/reference | Verdict |
|---|---|---|---|---|
| `apollo12_erosion_rate` | 0–16.6 kg/s (zero above its 31.9 m onset) | 10.7–99.0 kg/s | **0.00–0.43** | **fail** (band 0.33–3) |
| `apollo12_scour_radius` | 0–14.3 m (zero above its 31.9 m onset) | 5–68 m | **0.00–0.83** | **fail** (band 0.33–3) |
| `apollo_total_eroded_mass` | 622 kg | 2600 kg (1200–26,000 kg published) | 0.24 | pass (band 0.1–10) |
| `apollo_lm_surface_shear` | 0.17–6.11 Pa peak | 0.155–3.36 Pa area-average | 0.69–1.82 | pass (band 0.1–10) |
| `pfgt1_erosion_threshold_shear` | 0.150 Pa | 0.25 Pa (0.47 Pa on BP-1) | 0.60 | pass (band 0.2–5) |
| `surveyor3_ejecta_speed` | 140 m/s at 5 m | 300 m/s (range 40–2000) | 0.47 | pass (band 0.13–6.7) |
| `apollo16_erosion_onset_height` | 33.5 m | 31.5 m | 1.06 | calibration target — never gated |
| `viking_deep_cratering_pressure` | 5321 Pa at 1.83 m | 300 Pa limit | 17.7 | cross-regime — never gates |
| `apollo_dust_ejection_angle` | no output | 1.4–8.1°, mean 3.3° | — | not modeled |
| `apollo_scour_depth` | no output | 6 cm (3–12 cm) | — | not modeled |
| `ground_effect_augmentation` | 3.0% of thrust at contact | none published | — | unsourced |

## The integrated model

The model now integrates the erosion regimes of `regolith_erosion.jl` over the
footprint on the gas state the plume field reports at every radius, keeps a
radial grid of eroded depth under the vehicle, and evaluates the grain transport
of `ejecta_transport.jl` once per saved sample. The two fitted constants are
gone from the default path: the viscous regime takes its threshold from Shao and
Lu's expression applied to the Lunar Sourcebook's soil (0.057 Pa, against the
fitted 0.15 Pa), and no coefficient in it is tuned to anything in this study.

Two of the eleven cases changed status as a result. `apollo_scour_depth` became
gate-eligible, because the model now produces a crater depth. `apollo_dust_ejection_angle`
went from `not_modeled` to `calibration_target`, because the model now produces
an ejection angle but it is an input taken from that very table.
`apollo16_erosion_onset_height` went from `calibration_target` to
`different_quantity`: nothing is fitted to it any more, but the reference is the
height at which dust becomes *visible* and the model reports the height at which
soil starts to move, and those cannot be equal.

Both runs are the shipped defaults; nothing was tuned between them.

```bash
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --field=analytic
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --field=table
```

### Analytic field (the default): six of seven gate-eligible cases pass

| Case | Model | Reference | Model/reference | Verdict |
|---|---|---|---|---|
| `apollo12_erosion_rate` | 8.9–29.7 kg/s, non-zero at every measured altitude | 10.7–99.0 kg/s | 0.30–0.97 | **fail** (band 0.33–3, on the 1.83 m row alone) |
| `apollo12_scour_radius` | 2.5–23.0 m | 5–68 m | 0.34–0.91 | **pass** (was fail) |
| `apollo_total_eroded_mass` | 1200 kg | 2600 kg (1200–26,000 published) | 0.46 | pass |
| `apollo_lm_surface_shear` | 0.17–6.11 Pa peak | 0.155–3.36 Pa | 0.69–1.82 | pass |
| `apollo_scour_depth` | 3.19 cm deepest, at r = 0.18 m | 6 cm (3–12 cm) at 1–2 m | 0.53 | **pass** (newly gate-eligible) |
| `pfgt1_erosion_threshold_shear` | 0.0571 Pa | 0.25 Pa (0.47 Pa on BP-1) | 0.23 | pass (band 0.2–5) |
| `surveyor3_ejecta_speed` | 140 m/s at 5 m | 300 m/s (range 40–2000) | 0.47 | pass |
| `apollo16_erosion_onset_height` | 54.3 m | 31.5 m visible onset | 1.72 | different quantity — never gates |
| `apollo_dust_ejection_angle` | 1.95° | 2.6° (1.4–8.1°, mean 3.3°) | 0.75 | calibration target — never gates |
| `viking_deep_cratering_pressure` | 5321 Pa at 1.83 m | 300 Pa limit | 17.7 | cross-regime — never gates |
| `ground_effect_augmentation` | 3.0% of thrust at contact | none published | — | unsourced |

### Tabulated field (`data/psi/apollo_lmde.json`): five of seven pass

| Case | Model | Reference | Model/reference | Verdict |
|---|---|---|---|---|
| `apollo12_erosion_rate` | 24.1–28.5 kg/s | 10.7–99.0 kg/s | 0.28–2.26 | **fail** (band 0.33–3, on the 1.83 m row alone) |
| `apollo12_scour_radius` | 6.6–39.6 m | 5–68 m | 0.57–1.76 | pass |
| `apollo_total_eroded_mass` | 1396 kg | 2600 kg | 0.54 | pass |
| `apollo_lm_surface_shear` | 1.11–37.7 Pa peak | 0.155–3.36 Pa | 4.5–11.2 | **fail** (band 0.1–10) |
| `apollo_scour_depth` | 3.40 cm deepest, at r = 0.22 m | 6 cm (3–12 cm) | 0.57 | pass |
| `pfgt1_erosion_threshold_shear` | 0.0571 Pa | 0.25 Pa | 0.23 | pass |
| `surveyor3_ejecta_speed` | 75.7 m/s at 5 m | 300 m/s (range 40–2000) | 0.25 | pass |
| `apollo16_erosion_onset_height` | 140 m | 31.5 m visible onset | 4.45 | different quantity — never gates |
| `apollo_dust_ejection_angle` | 1.86° | 2.6° | 0.72 | calibration target — never gates |
| `viking_deep_cratering_pressure` | 6485 Pa at 1.83 m | 300 Pa limit | 21.6 | cross-regime — never gates |
| `ground_effect_augmentation` | 3.0% of thrust at contact | none published | — | unsourced |

## The calibration question the field work left open

The plume-field work reported that swapping the Gaussian footprint for the
tabulated source-flow plume moved the Apollo 11 descent from 429 kg of eroded
soil to 29 tonnes — a factor of 68 — because `threshold_shear_pa` and
`erosion_efficiency` had been calibrated against the analytic field's much
weaker wall shear and did not transfer.

**That problem is gone, and it needed no recalibration.** The two constants are
not used by the default law at all. The threshold is now a property of the soil
(`shields_threshold_shear_pa`, 0.0571 Pa on lunar mare at 1.625 m/s²) and does
not change when the field changes, so the same soil erodes under both plumes.
Total eroded mass over Lane and Metzger's Apollo 12 descent profile:

| | analytic | table | ratio |
|---|---|---|---|
| fitted Roberts law (the old default) | 622 kg | 16.8 t | 26.9 |
| erosion regimes (the current default) | 1200 kg | 1396 kg | **1.16** |

(The factor of 68 the field work quoted is the same effect measured over the
Apollo 11 descent the demo flies rather than over Lane and Metzger's Apollo 12
profile; both integrations show the same thing.)

The one constant that remains unsourced is the soil's `saltation_efficiency`,
the factor of 10 that stands in for the cascade of grains splashed up by each
impacting grain. It was **not** recalibrated here, and the reasoning is on the
record:

- The best fit to Lane and Metzger's eleven Apollo 12 altitudes is 17.5 on the
  analytic field and 13.7 on the table, both within a factor of 1.8 of the
  shipped 10. The shipped value therefore sits at 0.57 and 0.73 of the best fit.
- Fitting it to `apollo12_erosion_rate` would make that case circular, and this
  study's own rule 2 forbids gating on a calibration target. Refitting would
  convert the study's strongest gate into a tautology in exchange for a factor
  of 1.8.
- It would not fix the shape in any case. After rescaling by the best fit the
  residuals still run 0.53 to 1.69 on the analytic field and 0.38 to 3.11 on
  the table: a single multiplier cannot move a rate that is nearly flat with
  height onto a reference that rises 9.3-fold toward the ground.

So it is left at 10 and reported as the model's remaining fitted constant.

## What the model gets right and wrong

Written from the two runs above, not from expectation.

### Right

- **It erodes soil everywhere the measurement does.** The single biggest defect
  of the fitted law was that it returned exactly zero above 31.9 m while Lane
  and Metzger measure 10.7 kg/s at 36.6 m. The derived threshold puts the onset
  at 54 m (analytic) and 140 m (table), so every one of the eleven measured
  altitudes now has a non-zero model rate, and the worst ratio is 0.30 instead
  of 0.00.
- **The eroding region is the right width.** `apollo12_scour_radius` was the
  other failing case and it now passes on both fields: 0.34–0.91 on the analytic
  field, 0.57–1.76 on the table. The table is the better of the two, because its
  wall jet spreads shear far outside the momentum core, which is the geometry
  Lane and Metzger measure (an effective half-angle of 52° to 66°).
- **The two fields agree on how much soil moves.** 1200 kg against 1396 kg,
  where the old law differed by a factor of 68. The plume model and the erosion
  model are genuinely separable now.
- **The crater depth lands inside the only published estimate.** 3.2 cm
  (analytic) and 3.4 cm (table) against Metzger's 6 cm, inside his stated 3–12 cm
  extremes of plausibility. This is a new output, so it is a prediction rather
  than a fit — but it is at the shallow edge of the band, and the case passes at
  0.53 against a tolerance floor of 0.50.
- **Total eroded mass improved from 0.24 to 0.46 of the reference** without
  anything being tuned, simply from eroding over a wider region for longer.

### Wrong

- **The rate is too flat with height, and the table is worse.** The measured
  rate rises by a factor of 9.3 from 10.7 kg/s at 36.6 m to 99.0 kg/s at 1.83 m.
  The analytic field rises by 3.35 over the same span; the table is essentially
  flat (24.1 to 27.7 kg/s, a factor of 1.15, which is why its worst ratios sit at
  both ends of the profile instead of one). Both fail the lowest altitude, where
  the reference rises 63 percent in the last four meters and neither model rises
  at all.
- **The tabulated field's wall shear is four to eleven times the Apollo CFD
  fit.** `apollo_lm_surface_shear` passed on the analytic field and fails on the
  table at every height. Roberts' rough-surface drag coefficient of 0.2 acting on
  the wall-jet dynamic pressure is a far stronger coupling than a 0.01
  skin-friction coefficient on the static pressure, and Metzger, Lane, Immer and
  Li (Earth & Space 2010) say directly that the smooth-wall CFD the reference
  fits understates the real shear while Roberts' rough-wall approximation is
  orders of magnitude higher. The reference and the table bracket the truth
  rather than one of them being right, and this study cannot say where in
  between it falls.
- **The derived threshold is further from NASA's measurement than the fitted one
  was.** `pfgt1_erosion_threshold_shear` went from 0.60 to 0.23 of the 0.25 Pa
  Stubbs and Mehta read at the crater edge, because Shao and Lu's cohesion
  parameter, fitted to terrestrial soils, gives 0.057 Pa for lunar fines.
  Matching 0.25 Pa would need that parameter 2.8 times the top of their fitted
  range — plausible for airless, electrostatically charged fines, but unmeasured,
  so it is left alone. The case still passes its factor-of-five band.
- **The erosion onset is high, and on the table it is very high.** 54 m and
  140 m against a 31.5 m *visible* onset. The physical onset must be above
  36.6 m, where erosion is still measured, so 54 m is not obviously wrong; 140 m
  is hard to defend and follows directly from the table's shear being four to
  eleven times the CFD fit.
- **The crater depth is an upper bound, which makes its agreement worse than it
  looks.** The model credits every second of erosion to one patch of ground: the
  grid is axisymmetric about the current impingement point and does not follow
  it. On the Apollo 11 descent the vehicle covers 310 m horizontally below the
  erosion onset height, and 63 percent of the eroded mass falls in the last
  92 m of it, against a crater edge radius of 4.3 m. So 3.2 cm is what a
  *hovering* vehicle would dig over that erosion history, and the real
  localized depth is lower — the case passes at 0.53 of a reference it is
  already under-predicting, and the true residual is larger than that.
- **The crater's deepest point is in the wrong place.** The model puts it at
  0.18–0.22 m from the centerline; Metzger puts it at 1–2 m. The depth is right
  to a factor of two and the shape is toroidal as reported, but the torus is
  five to ten times too tight. This is the same "the plume is concentrated into
  too small a patch" defect the baseline had, seen in a new quantity, and it is
  why the ejection angle is *not* wired to the crater's wall slope: the slope
  this profile gives is far steeper than the 1–3° measured.
- **The deposition radius is hundreds of meters to kilometers,** against tens of
  meters of visible blast zone at the Apollo sites. The fix the ejecta work
  identified — weight the launch radii by the regimes' own local erosion rate
  instead of by the raw wall shear, and weight the grain sizes by a lognormal
  mass distribution fitted to the soil's D50 and D84/D50 — halves it and no
  more. At the Apollo approach thrust 10 m up, the mass-weighted mean deposition
  radius goes from 1553 m to 755 m and the mean speed from 143 to 107 m/s; the
  escape fraction is 0.0000 either way on this plume field, at both 11.5 kN and
  45 kN. The remainder is not a weighting artifact: a 70 µm grain leaving at
  107 m/s and 2° above the horizontal has a ballistic range of 424 m in lunar
  gravity, and that follows directly from the measured ejection angle and the
  measured escape speeds. The visible blast zone is where the coarse mass lands,
  not where the population lands; the model has no separate output for it and
  the two should not be compared.
- **The ground-effect force still has no reference at all.** Unchanged: three
  percent of thrust at contact decaying over two exit diameters is an
  assumption, the literature has not characterized the near-ground base-pressure
  rise, and this study cannot say whether it is right to a factor of ten. It is
  also the only wrench the effector applies to the vehicle.

### The shape of the gap

The baseline's diagnosis was that one geometric error — a plume spread over
`h·tan 25°` where the soil moves out to `1.3 h` to `2.2 h` — was hiding inside a
calibrated constant. That is fixed: the width now passes on both fields and the
calibrated constants are out of the default path. What is left is a different
error, and it is in the *height dependence*. Three published laws disagree about
it and this model agrees with none of them:

| Source | Dependence of the surface stress on height |
|---|---|
| this model, both fields | `h^-2` (a point source on a plane) |
| Morris 2012, DSMC of the LMDE | `h^-2.825` |
| Lane and Metzger 2015, Apollo 12 film fit, Eq. (14) | `exp(-0.123 h)` |

Over 1.83 m to 36.6 m those three span a factor of 30 in how much the stress
grows as the vehicle descends, and the erosion rate follows the stress. Nothing
in this study can choose between them: the Apollo film fit is an exponential
with no physical derivation, the DSMC is a single engine at a single condition,
and `h^-2` is what a point source must give. It is recorded as the open question
it is.
