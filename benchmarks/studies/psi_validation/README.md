# Plume-surface interaction validation

This study evaluates SpaceAGORA's plume-surface interaction model
(`src/dynamics/coupled/force_torque_models/plume_surface_interaction.jl`,
documented in [`docs/src/user/lunar_landing.md`](../../../docs/src/user/lunar_landing.md))
against every published plume-surface measurement that could be sourced
precisely, and records the gaps.

It is modeled on [`benchmarks/studies/telemetry_validation/`](../telemetry_validation/README.md):
one frozen manifest per reference case under `manifests/`, a runner that writes
a comparison table to `results/` (gitignored), and `common.jl` for the shared
parts. Unlike that study it runs no simulation — every case is a closed-form
evaluation of the model's public functions (`plume_surface_footprint`,
`plume_quantities`, `plume_erosion_onset_height`, `plume_ground_effect_force`),
so the whole thing takes under a second after Julia loads, needs no SPICE, no
GRAM and no results file.

## Commands

```bash
# Every case, printed and written to results/<hostname>/ (nothing gates)
julia --project=. benchmarks/studies/psi_validation/run_validation.jl

# Only the cases that may gate, and fail the run if one is outside tolerance
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --cases=gates --enforce=true

# One case
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --cases=pfgt1_erosion_threshold_shear
```

Outputs (gitignored):

```
results/<hostname>/
  psi_validation_comparison.csv   # one row per reference value: model, reference, ratio, verdict
  psi_validation_cases.csv        # one row per case: status, gating, tolerance, full citation
```

The run ends with a findings block computed from the rows it just produced: the
model/reference span of every case, whether the model's trend with height has
the same sign as the reference's, where each series peaks, and what is on record
as unmeasurable. The prose below is written from that block and the table above
it, not from expectation.

Three of the cases also run as PR-tier regression tests in
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
   loader refuses any other combination. The four non-gating statuses are
   `cross_regime` (a real measurement from a different body or flow regime),
   `calibration_target` (the observable a model parameter was fitted to, so
   agreement is circular), `not_modeled` (the model produces nothing comparable)
   and `unsourced` (no published value was found at all).
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
| `apollo16_erosion_onset_height` | height at which dust starts blowing | Metzger 2024 part 2, §2 | no — calibration target |
| `viking_deep_cratering_pressure` | surface pressure at which deep cratering starts | Metzger 2016, "Deep cratering on Mars" | no — cross-regime |
| `apollo_dust_ejection_angle` | ejection angle above horizontal | Immer et al. 2008, Table 1 | no — not modeled |
| `apollo_scour_depth` | depth of the toroidal scour | Metzger 2024 part 2, §5 | no — not modeled |
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
  case. The model has no crater depth output in any event.
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

## Baseline: today's model

Recorded on `space-falcon-1` from commit `5c4ed090`'s model, with
`PlumeSurfaceConfig()` defaults (Apollo LM DPS, 1.5 m exit diameter, 25°
half-angle, `friction_coefficient` 0.01, `threshold_shear_pa` 0.15,
`erosion_efficiency` 10, lunar mare soil). **Four of six gate-eligible cases
pass; two fail.** These numbers are the reference point the next round's model
has to beat.

| Case | Model | Reference | Model/reference | Verdict |
|---|---|---|---|---|
| `apollo12_erosion_rate` | 0–16.6 kg/s (zero above its 31.9 m onset) | 10.7–99.0 kg/s | **0.00–0.43** | **fail** (band 0.33–3) |
| `apollo12_scour_radius` | 0–14.3 m (zero above its 31.9 m onset) | 5–68 m | **0.00–0.83** | **fail** (band 0.33–3) |
| `apollo_total_eroded_mass` | 622 kg | 2600 kg (1200–26,000 kg published) | 0.24 (0.024–0.057 vs. the 2024 estimate) | pass (band 0.1–10) |
| `apollo_lm_surface_shear` | 0.17–6.11 Pa peak | 0.155–3.36 Pa area-average | 0.69–1.82 | pass (band 0.1–10) |
| `pfgt1_erosion_threshold_shear` | 0.150 Pa | 0.25 Pa (0.47 Pa on BP-1) | 0.60 (0.32) | pass (band 0.2–5) |
| `surveyor3_ejecta_speed` | 140 m/s at 5 m | 300 m/s (range 40–2000) | 0.47 | pass (band 0.13–6.7) |
| `apollo16_erosion_onset_height` | 33.5 m | 31.5 m | 1.06 | calibration target — never gates |
| `viking_deep_cratering_pressure` | 5321 Pa at 1.83 m | 300 Pa limit | 17.7 | cross-regime — never gates |
| `apollo_dust_ejection_angle` | no output | 1.4–8.1°, mean 3.3° | — | not modeled |
| `apollo_scour_depth` | no output | 6 cm (3–12 cm) | — | not modeled |
| `ground_effect_augmentation` | 3.0% of thrust at contact | none published | — | unsourced |

## What the model gets right and wrong

Written from the run above, not from expectation.

### Right

- **The shear stress it puts on the ground is the right size.** The peak wall
  shear runs 0.69 to 1.82 times the Apollo LM CFD fit across 5 to 30 m — the
  best agreement in the study, and over a factor of twenty in the reference
  value. Given that the model reaches that number through a guessed Gaussian
  footprint and a single skin-friction coefficient, landing inside a factor of
  two of a RANS CFD is more than it had any right to.
- **The erosion threshold is defensible.** The 0.15 Pa fitted constant, chosen
  to place the onset at 31 m, turns out to sit at 0.60 of the 0.25 Pa NASA
  measured at the crater edge in a vacuum chamber, and 0.32 of the 0.47 Pa on
  the coarser lunar simulant. A number fitted to a 1969 crew comment landed
  within a factor of two of a 2026 measurement of the same physical threshold.
  This was not known when the constant was chosen and is the strongest single
  result here.
- **The characteristic ejecta speed is inside the measured range.** 140 m/s at
  5 m sits between the Surveyor III lower bounds (40, 70, 100 m/s) and the
  pitting-morphology estimate (300–2000 m/s). It is at the low end, and the
  compilation calls the pitting estimate the most reliable, so "inside the
  range" is all this case can claim.
- **The onset height matches the measurement** — but that case is a calibration
  target and proves nothing. It is recorded so the circularity stays visible.

### Wrong

- **It moves far too little soil.** Against
  Lane and Metzger's Apollo 12 profile the model returns 0.10 to 0.43 of the
  measured rate at every altitude below its onset, and nothing at all above it;
  the deficit is worst nearest the ground (0.14 at 1.83 m) and just under the
  onset (0.10 at 29.3 m). The total over the descent is 622 kg against their 2600 kg (0.24), or 0.024 to
  0.057 of the 11–26 t Metzger published in 2024. The lumped
  `erosion_efficiency = 10` is absorbing a factor that is really four to forty.
- **The rate has its maximum in the wrong place.** The measured rate rises as
  the vehicle descends, from about 11 kg/s at 36.6 m to 99 kg/s at 1.83 m, with
  row-to-row scatter of order 40 percent on the way but its maximum at the
  lowest altitude. The model's rate instead peaks at about 16.6 kg/s near 11 m and then
  *falls* to 14.2 kg/s at contact. The cause is structural: the Gaussian
  footprint is normalized so that `∫p dA = F` at every height, so as the
  footprint tightens the excess-shear integral saturates and then shrinks with
  the area. No amount of retuning the two constants fixes the sign of that
  slope.
- **The eroding region is much too narrow, and increasingly so with height.**
  The measured `a0/h` stays between 1.3 and 2.2 over most of the descent (2.7 at
  the lowest altitude sampled) — an effective half-angle of 52° to 66°, far wider than the
  25° momentum core the model spreads pressure over. The model's annulus is
  1.27 h at 1.83 m but only 0.65 h at 21 m and 0.44 h at 30 m, because the outer
  shear falls through the threshold; and it could not exceed 1.4 h even if the
  shear held up, because its quadrature stops at three footprint radii. Against
  the reference the ratio therefore decays from 0.83 to 0.25.
  The wall jet spreads soil far outside the momentum-carrying core, and this
  model has no wall jet.
- **Erosion switches off abruptly at 31.9 m, where the measurement is still
  reading 10.7 kg/s.** The hard threshold gives exactly zero above the onset
  height, and the highest altitude Lane and Metzger analyzed is above it. The
  real onset is higher still: they see erosion at the top of their usable range,
  not at its edge.
- **The area-averaged shear is a factor of four low, which means the peak
  agreement is partly luck.** Averaging the model's own shear profile over the
  reference's erosion radius the way the reference defines its average gives
  0.42 Pa at 10 m against 1.82 Pa (0.23), and 0.097 Pa at 15 m against 0.98 Pa
  (0.10). The peak matches because the model concentrates its shear in a narrow
  ring; spread over the area where soil actually moves, there is four to ten
  times too little of it. That is the same defect as the narrow footprint, seen
  in a different quantity, and it is consistent with the erosion-rate deficit.
- **The surface pressure exceeds the only published deep-cratering limit by a
  factor of eighteen, at a landing that dug no crater.** 5.3 kPa at 1.83 m
  against the Viking 0.3 kPa empirical limit. The limit is Martian and its own
  source warns against transferring it, so this case cannot gate — but it points
  the same way as everything else: the model concentrates the plume's momentum
  into too small a patch.
- **The ground-effect force — the only wrench the effector actually applies to
  the vehicle — has no reference at all.** Three percent of thrust at contact
  decaying over two exit diameters is an assumption, the literature has not
  characterized the near-ground base pressure rise, and this study cannot say
  whether it is right to a factor of ten.

### The shape of the gap

Five of the seven defects above are the same defect: the plume is spread over
`h·tan 25°` when the soil moves out to `1.3 h` to `2.2 h`. Concentrating the
thrust into a patch four to twenty times too small in area raises the peak
pressure and peak shear into agreement while starving the area integral that
sets the erosion rate — which is then partly compensated by the
`erosion_efficiency = 10` multiplier, hiding the geometric error inside a
calibrated constant. A plume field computed from the nozzle with a real wall
jet, rather than a normalized Gaussian, is what these numbers ask for.
