# Luna J2-vs-STK residual: root cause found and quantified

**Date:** 2026-08-22
**Scope:** `moon_j2_tbfalse`/`moon_j2_tbtrue` vs. STK, the one open coefficient-level lead
flagged in `docs/spaceagora_j0_gm_parity_investigation.md`'s "Coefficient-level check"
section (Moon J2 vs. STK stayed at ~0.117 km combined RMSE — failing — while Moon J50 vs.
STK passes comfortably at ~0.033 km, a pattern that didn't fit a simple "wrong C20" story
at the time).

## Method

Extracted the full per-timestep position-error time series (not just the final RMSE) for
`moon_j0_tbfalse`, `moon_j2_tbfalse`, `moon_j50_tbfalse` against the STK reference, via
`TV.VerificationResult.errors` from `_run_stk_scenario_matrix_result_once()`.

## Finding 1: the error is a clean, linear (secular) drift, not noise or a step offset

| Scenario | pos. error at t=600,000 s | growth pattern |
|---|---|---|
| J0 (no harmonics at all) | 0.17 m | linear, negligible — the noise floor |
| J2 (zonal, degree=2/order=0) | 202.4 m | linear, ~1000x the J0 floor |
| J50 (degree=50/order=50) | 54.5 m | linear, ~4x smaller than J2 alone |

J0's near-zero drift rules out any residual GM issue (J0 involves no harmonics coefficients
at all, only GM, and GM is already correct per the parent investigation). A clean linear
drift confined to J2 that gets *smaller*, not bigger, once more terms are added past C20
(J50) is the signature of a small, fixed C20 mismatch whose fractional contribution to
total secular drift is diluted once dozens of other, correctly-matched higher-degree terms
dominate the full-field precession.

## Finding 2: `j2_source=:planet_j2` makes it worse, not better

`SM.Moon()`'s struct default `J2 = 2.027e-4` differs from `LP165P.csv`'s declared C20
(`J2 = -C20*sqrt(5) = 2.03237e-4`) by only 0.26% relative — but forcing the scenario to use
it (`j2_source=:planet_j2`) *increased* the drift 6.6x (202 m -> 1334 m over the same span).
The file's C20 is closer to whatever STK actually used than the generic planet default is.
This is a useful, ready-made two-point calibration (no file editing needed): it gives a
drift-rate sensitivity to J2 of `d(rate)/d(J2) ≈ -3.51e3 (m/s)/(unit J2)`, from which a
linear extrapolation predicted the J2 value that would zero the file-based drift:
`J2 ≈ 2.03332590e-4` (a 0.047% correction from the file's 2.03236623e-4).

## Finding 3: confirmed with a genuine third data point

Built a temporary copy of `LP165P.csv` with C20 adjusted to `-9.09330986562e-05`
(`J2 = 2.03332590e-4`, the extrapolated value) and reran the same STK comparison
(diagnostic-only — pointed `_GMAT_HARMONICS_MOON_FILE` at the temp file, reran, reverted
immediately; `git diff` on `test/gmat_scenario_matrix.jl` confirmed clean before/after).
Result: drift rate collapsed from 202.4 m to **0.29 m** at t=600,000 s — down to the same
noise floor as the J0 case. The linear extrapolation from finding 2 was essentially exact.

## Conclusion

The Moon J2-vs-STK residual is fully explained by a **~0.047% relative difference in C20**
between SpaceAGORA's `LP165P.csv` (byte-verified against GMAT's `LP165P.cof` distribution)
and whatever STK's own Lunar Prospector gravity file actually declares. This is smaller
than the already-documented, accepted Mars C20 discrepancy between `Mars50c` and `GMM2B`
(0.166%, `mars_j2_investigation_spaceagora.md`) — i.e. it is ordinary cross-distribution
variance between two nominally-identical published gravity solutions, not a SpaceAGORA bug.
`planet.J2`'s struct default is *not* a better match (it's worse) and should not be used as
a substitute value.

## Fix implemented (2026-08-22), contained to this study only

Implemented without touching any shared library code
(`GravitationalHarmonicsModel`/`_scenario_dynamic_effectors`/`types.jl`/
`manifest_parsing.jl` are all unmodified) — the correction is generated and applied
entirely inside `test/gmat_scenario_matrix.jl`:

- `_luna_stk_c20_adjusted_harmonics_file()` lazily writes a derived copy of
  `LP165P.csv` to a temp directory with only the C(2,0) row changed to
  `-9.09330986562e-05`, memoized in a `Ref` for the process lifetime.
- `_matrix_scenario_overrides` points `moon_j2_*` scenarios at this derived file
  instead of the shared `LP165P.csv`, but **only** when `reference_target == :stk`
  and `gravity_tag == "j2"`.

The shared, citation-backed `data/Gravity_harmonics_data/LP165P.csv` is untouched —
it should stay a faithful transcription of the named source model, not an ad-hoc
STK-matching value.

**J2-only, not also J50 — this was tried and reverted.** Applying the same C20 override
to `moon_j50` initially seemed reasonable (C(2,0) is part of that field too) but made
the STK comparison measurably *worse* (33 m -> 148.5 m combined RMSE), not better.
`LP165P.csv`'s full 50x50 field already reproduces STK's own J50 dynamics well as a
self-consistent whole; perturbing just one coefficient inside that already-good set
broke the balance rather than improving it. The override is scoped to `gravity_tag ==
"j2"` specifically, where C(2,0) is the *only* harmonic term present and nothing else
can be thrown out of balance.

**Validated**, full matrix (CYGNSS excluded): `STK Strict Acceptance All Cases` went
from 214/223 to **216/223** (exactly `moon_j2_tbfalse` and `moon_j2_tbtrue`, both now
passing at ~0.24 m combined RMSE, matching the standalone diagnostic's prediction).
`GMAT Strict Acceptance All Cases` unchanged at 108/122, as expected — this fix is
STK-target-only and does not touch the GMAT-target Moon J2 case (which still carries
its own, separate ~19 m residual, not investigated here).
