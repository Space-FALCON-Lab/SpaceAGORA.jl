# J0/J2/J50 GM and gravity-order parity: SpaceAGORA vs. GMAT vs. STK

**Date:** 2026-08-20
**Scope:** all `<body>_<j0|j2|j50>_<tbfalse|tbtrue>` scenarios in the GMAT/STK scenario
matrix (`test/gmat_scenario_matrix.jl`). CYGNSS scenarios in the same file are out of
scope for this document.

## Question

For the J0 (pure point-mass, degree/order 0) case of each body, is the central-body
gravitational parameter (GM) SpaceAGORA uses for its own simulated trajectory the same
value the checked-in GMAT and STK reference trajectories were actually generated with?

## Method

Two-body (Keplerian) motion conserves specific angular momentum `h = r × v` and the
semi-latus rectum `p = a(1-e²) = h²/μ` exactly. Both are queryable straight from a single
row of a reference trajectory (`r`, `v`, and the file's own reported `SMA`/`ECC` columns),
so `μ = h² / (a(1-e²))` recovers the *exact* GM the propagator that generated the reference
data used — no full-trajectory simulation or curve-fitting required, and (unlike
differencing radii across samples) it is numerically well-conditioned even for the
near-circular J0 orbits used here (e as low as 0.00001).

Applied to:
- `data/telemetry/Basilisk_Examples_Full/Sim_{Body}_1M_J0_TBFalse.feather` — the file
  `test/gmat_scenario_matrix.jl`'s `_GMAT_EXAMPLES_DIR`/`_scenario_basilisk_path` actually
  resolves to for the "GMAT Strict Acceptance" testset (despite the directory's name).
- `data/telemetry/stk_results/Sim_{Body}_J0_TB0.feather` — the STK reference.

## Finding: GMAT's and STK's own J0 references don't agree with each other

Reconstructed GM (km³/s²), 5 checks per body, all agreeing to displayed precision:

| Body  | GMAT-reference μ | STK-reference μ | Relative difference |
|-------|-----------------:|-----------------:|---------------------|
| Earth | 398600.436000     | 398600.441500     | 1.4e-8 |
| Mars  | 42828.314258      | 42828.372854      | 1.4e-6 |
| Venus | 324858.599000     | 324858.589726     | 2.9e-8 |
| Moon  | 4902.799000       | 4902.800306       | 2.7e-7 |

The GMAT-reference values for Earth/Venus/Moon are exactly SpaceAGORA's own bundled
`data/GRAMSuite.jl/GRAM Suite 2.0/SPICE/pck/de_403_masses.tpc` kernel constants
(`BODY399_GM=398600.436`, `BODY299_GM=324858.599`, `BODY301_GM=4902.799`) — i.e.
SpaceAGORA's existing `planet.μ` already matched the GMAT reference almost exactly for
those three bodies. **Mars was the sole outlier**: SpaceAGORA forces
`MARS_MU_M3S2 = 42828.372854` (`src/environment/ephemerides/planets.jl`), a value that
matches the *STK* reference almost to machine precision but differs from the GMAT
reference by ~0.06 km³/s² — small in relative terms, but enough to dominate the J0
position-error budget over a 1,000,000 s / ~150-orbit propagation.

Confirmed empirically by temporarily swapping Mars's `μ` to the GMAT-reference value and
re-running the matrix: the GMAT-comparison RMSE dropped 16× (0.374 km → 0.023 km) while
the STK-comparison RMSE *grew* from ~5e-6 km (near machine precision) to ~2.2 km. There is
no single Mars `μ` that satisfies both references simultaneously — this is not a
SpaceAGORA bug so much as GMAT's and STK's own J0 scenario generators using different GM
conventions for the same body (GMAT's generator skips loading a potential file at
degree=0 and falls back to its own default per-body Mu; STK's generator always loads the
gravity file, even at degree=0, and just truncates the expansion).

## Fix

Rather than pick one convention and accept a mismatch against the other tool, each
comparison target now runs SpaceAGORA's J0 case with the density-matched GM for *that*
target specifically:

- `TimeAlignedScenarioConfig`/`OrbitEventsScenarioConfig` gained an optional
  `gravity_harmonics_gm_override_m3s2` field
  (`src/analysis/verification/telemetry_verification/types.jl`, wired through
  `manifest_parsing.jl`).
- `_scenario_dynamic_effectors` (`scenario_builders.jl`) builds a degree-0
  `GravitationalHarmonicsModel` with that GM explicitly passed as `gm_m3s2` whenever the
  override is set, instead of falling through to `InverseSquaredGravityModel` (which reads
  the generic `planet.μ`). When the override is absent, behavior is unchanged.
- `test/gmat_scenario_matrix.jl`'s `_matrix_scenario_overrides` takes a
  `reference_target::Symbol` (`:gmat` or `:stk`) and sets the override from the table above
  for `*_j0_*` scenarios; `_run_gmat_scenario_matrix_result_once`/
  `_run_stk_scenario_matrix_result_once` each run their own independent SpaceAGORA
  simulation (they already did — separate manifests, separate `TV.run_verification` calls)
  with `reference_target=:gmat` / `:stk` respectively.

J2/J50 scenarios are untouched: their GM already comes from each harmonics file's own
declared `gm_m3s2` header, not `planet.μ`, and those files are already byte-identical to
the GMAT `.cof` distributions for Mars/Venus/Moon (verified header-for-header).

## Result

J0 strict-acceptance RMSE, `tbfalse`, before → after (km):

| Body  | vs. GMAT | vs. STK |
|-------|----------|---------|
| Earth | 0.052 → 0.052 (unchanged — dominated by something else, not GM) | 0.046 → 5.9e-6 |
| Mars  | 0.389 → 0.023 | 5e-6 → 5.0e-6 (unchanged) |
| Venus | 8.3e-7 → 8.3e-7 (unchanged) | 0.095 → 6.3e-6 |
| Moon  | 1.5e-5 → 1.5e-5 (unchanged) | 0.205 → 6.8e-5 |

Full matrix (`GMAT Strict Acceptance All Cases` + `STK Strict Acceptance All Cases`,
all bodies × J0/J2/J50 × tbfalse/tbtrue) run clean afterward: STK 214/223 passing, GMAT
104/122 passing. Remaining failures are in J2/J50/`tbtrue` cases unrelated to this fix
(known pre-existing issues — see `mars_j2_investigation_spaceagora.md` for the J2 tesseral
finding) plus one CYGNSS legacy borderline case using degree=50, not touched by this
change.

Earth's and Mars's remaining ~0.02-0.05 km GMAT-comparison residual is *not* a GM issue
(GM is now dialed in to ~1e-8 relative precision) — it's caused by something else specific
to the GMAT/`Basilisk_Examples_Full` reference for those two bodies. Worth a separate
investigation; out of scope here.

## Extension: J2 and J50

Applying the same treatment to J2/J50 turned out to need two separate fixes, not one —
GM alone was not the whole story here.

### J2 tesseral order: GMAT's and STK's "J2" are different fields

`_matrix_scenario_overrides` was fixed on 2026-08-19/20 (see `mars_j2_investigation_spaceagora.md`)
to request zonal-only (degree=2, order=0) for every body's J2 case, matching what both
tools' *scenario generator scripts* say they intend. But the checked-in
`Basilisk_Examples_Full` reference data was never regenerated after that fix — for
Mars/Venus/Moon it is still a full tesseral (degree=2, order=2) field (Earth's reference
already was zonal). Meanwhile STK's checked-in reference genuinely is zonal-only for every
body (confirmed via `generate_stk_cases.py`'s `GRAVITY_ORDERS = {0: 0, 2: 0, 50: 50}`).
So "GMAT J2" and "STK J2" are, as checked into this repo, two different physical models
for Mars/Venus/Moon.

Confirmed by temporarily forcing SpaceAGORA back to tesseral (2,2) for all bodies and
re-running just the J2 cases (`tbfalse`, combined-axis RMSE, km):

| Body  | zonal SA vs. GMAT ref | tesseral SA vs. GMAT ref | zonal SA vs. STK ref | tesseral SA vs. STK ref |
|-------|-----------------------|---------------------------|------------------------|----------------------------|
| Mars  | 50.6 | **0.54** | **2.2e-4** | 291 |
| Moon  | 8.5  | **0.07** | **8.3e-2** | 50.2 |
| Venus | 0.14 | **0.01** | **5.9e-6** | 0.88 |
| Earth | 0.05 | 0.55 (worse) | 0.003 | 3.4 (worse) |

Fix: `_matrix_scenario_overrides` now takes the J2 order from
`_MATRIX_J2_ORDER_OVERRIDE[(planet, reference_target)]` — 2 (tesseral) for
Mars/Venus/Moon against GMAT, 0 (zonal) against STK and for Earth against either.
Degree stays 2 for both targets.

### GM: GMAT uses one uniform per-body value across J0/J2/J50; STK does not

The J0 GM table (above) turned out to generalize cleanly to GMAT but not to STK:

- Applying the J0-reconstructed GM at J2/J50 for the **GMAT** target improved every case,
  confirming GMAT's reference data uses the same fixed per-body GM regardless of degree
  (Mars J50 0.274->0.085 km, Moon J50 0.068->0.016 km, Venus J50 0.0124->0.0016 km, on top
  of the J2 tesseral-order fix above).
- Applying the same J0-reconstructed GM at J2/J50 for the **STK** target made things
  *worse* (Venus J2 5.9e-6->0.024 km, Mars J2 2.2e-4->0.093 km) — STK's J2/J50 dynamics
  already match the harmonics file's own declared `gm_m3s2` almost exactly (that's why
  those cases were already excellent before any GM override existed), and only STK's J0
  case uses a distinct value.

Fix: the GM override applies at every degree for `:gmat`, but only at J0 for `:stk`
(`gravity_harmonics_gm_override_m3s2` is left unset — falls back to the file's own
`gm_m3s2` — for `:stk` J2/J50).

### Coefficient-level check (Clm/Slm, not just GM)

Given J2/J50 against STK were already excellent for Mars and Venus *before* any override
existed (Mars J50: 1.7e-4 km, Venus J50: 5e-6 km, both over ~150 orbits) and stayed
excellent after, that is itself strong evidence: any real Clm/Slm mismatch between
SpaceAGORA's/GMAT's `Mars50c.csv`/`MGNP180U.csv` and whatever STK's own `.grv` files
actually declare would show up as secular drift growing over the arc, and none is visible
at the sub-meter level for those two bodies at either degree.

Moon is the exception: `moon_j2_tbfalse` vs. STK stayed at ~0.117 km combined RMSE
throughout (unaffected by either fix above — it uses the same zonal order and file GM
before and after), while `moon_j50_tbfalse` vs. STK passes comfortably (0.033 km, under
the 0.05 km limit). That degree-2-specific-but-not-degree-50 pattern doesn't fit a simple
"STK's Moon C20 differs from LP165P.csv's" story cleanly (a wrong C20 would show up in the
J50 case too, since J50 includes C20). A secular-nodal-regression-based attempt to back out
STK's own apparent J2 from its Moon trajectory did show STK's implied J2 running ~4% away
from GMAT's implied J2 computed the same way (vs. ~2% for Mars) — consistent with *something*
being off specifically for the Moon — but that reconstruction has an unresolved reference-frame
caveat (computed inclination did not match the design value, e.g. 5.5 deg recovered for Mars
against a 45 deg design inclination, so the "PlanetInertial"/"Inertial" export frame's pole
is not what was assumed) and should be treated as suggestive, not confirmed. Root cause not
found; flagged here as the one open coefficient-level lead rather than chased further.

### Result after both fixes

Full matrix (`GMAT Strict Acceptance All Cases` + `STK Strict Acceptance All Cases`, all
bodies x J0/J2/J50 x tbfalse/tbtrue, CYGNSS excluded from this comparison):

| | before (J0 fix only) | after (J0+J2+J50 fix) |
|---|---|---|
| GMAT | 104/122 | **108/122** |
| STK  | 214/223 | **214/223** (no regression) |

Remaining GMAT failures: Earth J0/J2/J50 and Mars J2/J50 (the pre-existing, not-fully-explained
~0.02-0.3 km SA-vs-reference residual documented in `mars_j2_investigation_spaceagora.md`,
attributed there to differing source gravity coefficient files/reference radii rather than
gravity-model completeness), plus Moon J2. Remaining STK failures: `tbtrue` (third-body)
cases for Mars/Moon/Venus J0 — central-body GM/order is now correct, but third-body
(Sun/Earth) GM/ephemeris parity between SpaceAGORA and STK has not been checked and is a
plausible next lead — plus Moon J2 (see coefficient-level check above) and a narrowly-missed
Earth J50 case (0.059 km vs. 0.05 km limit).
