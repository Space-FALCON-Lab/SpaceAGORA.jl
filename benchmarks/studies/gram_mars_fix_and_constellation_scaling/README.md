# GRAM Mars Fix and Constellation Size Scaling

Two related pieces of work from the same investigation, kept together because
the scaling study depends on the Mars fix to be able to run Mars scenarios at
all.

## The Mars GRAM fix

`libGRAM.dylib` (the vendored native atmosphere library, `data/GRAMSuite.jl`)
statically links its own private copy of CSPICE, completely isolated from the
CSPICE instance SpaceAGORA's own `SPICE.jl` uses. GRAM's internal
`SpiceLoader` populates that isolated instance from its own default per-body
SPK kernels, which are Sun-inclusive only for Earth/Venus (a full planetary
ephemeris); Mars's default (`mar097_GRAM.bsp`) is Mars-system-only. As a
result, every Mars density query used to fail immediately with:

```
GRAM update failed (code=1): Error: A Spice error occurred.
       Error in longitude of the sun calculation.
```

This is not a missing-kernel-data problem -- the same `lspcn_c`/`ltime_c`
calls succeed when replayed against the identical kernel set via `SPICE.jl`.
It is that GRAM's own isolated CSPICE instance never gets fed the data.

The fix (`data/GRAMSuite.jl/src/GRAMSuite.jl`'s `_GRAM_EPHEMERIS_STATE_FN`
hook, populated by `ext/SpaceAGORAGRAMSuiteExt.jl`'s
`_gram_spice_ephemeris_state`) computes GRAM's required ephemeris state
(solar longitude, one-way light time, subsolar point, local solar time, solar
zenith angle, seconds-per-sol) using SpaceAGORA's own working `SPICE.jl`
bindings and feeds it to GRAM directly via `set_ephemeris_state!`, bypassing
GRAM's own broken internal SPICE calls entirely. The hook defaults to
`nothing` (GRAM's native behavior, unaffected) unless the extension is
loaded, and is registered for every body GRAM supports, not just Mars -- the
same isolated-CSPICE limitation plausibly affects Jupiter/Saturn/Uranus/
Neptune/Titan too (same satellites-only default kernel pattern), though only
Earth and Mars have been exercised end-to-end.

`mars_150km_gram_scaling.jl` is both a working example and a regression
check: if it starts failing again with the "longitude of the sun" error
above, the ephemeris bypass has broken.

```bash
julia --project=. --threads=4 benchmarks/studies/gram_mars_fix_and_constellation_scaling/mars_150km_gram_scaling.jl parallel
```

## Constellation size scaling

`leo_constellation_size_scaling.jl` sweeps LEO constellation size (1, 2, 4,
..., 1024 satellites -- powers of 2) at a fixed 600s mission, comparing:

- **standard**: real/native GRAM, no vacuum-predicted lookahead cache
  (`SPACEAGORA_DENSITY_FREEZE_PER_STEP=1`, density calls serialized)
- **surrogate**: the GRAM offline surrogate (`GRAMAtmosphereModelSurrogate`)
  -- precomputed interpolation grid, no native calls, no lock
- **no_gram**: `NoAtmosphereModel()` with the aerodynamic effector dropped
  entirely -- baseline gravity-only dynamics cost, isolating
  constellation-size scaling with zero atmosphere-model overhead

Each `(N_SATS, mode)` point runs in its own subprocess
(`leo_constellation_size_scaling_worker.jl`), rather than all 33 points in
one process, because `ODEParams` is parameterized on both `N_sats` and the
density-model type -- each distinct combination triggers its own JIT
specialization of the whole RHS/solver pipeline, and accumulating 33 of those
in one process risks the kind of memory pressure this project has repeatedly
hit on constrained machines.

```bash
julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_constellation_size_scaling.jl
```

Override the thread count used for every worker subprocess (default 4):

```bash
SPACEAGORA_SCALING_THREADS=8 julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_constellation_size_scaling.jl
```

Produces `leo_constellation_size_scaling_with_gram.png` (standard vs.
surrogate, both GRAM-derived) and `leo_constellation_size_scaling_without_gram.png`
(the no-GRAM baseline), each with a dashed O(N) linear-scaling reference line
anchored at N=1 for visually judging super-/sub-linear scaling.

**Key finding**: the lookahead/vacuum-predicted-cache GRAM mode -- used
elsewhere in this repo's benchmarks as the default "fast path" -- is actually
*slower* than the "standard" (no cache) mode at every constellation size
tested here (1 to 1024 satellites), because a 600s mission is too short to
amortize the cache's upfront drag-free-trajectory + spline-fit cost. The
lookahead cache's benefit is framed around much longer missions elsewhere in
this repo (e.g. a 1-hour mission) -- check actual mission duration before
assuming the lookahead cache is the fast option.

Companion single-point scripts used to establish this before building the
full sweep:

- `earth_standard_gram_scaling.jl` -- standard (no-lookahead) real GRAM, one
  constellation size at a time (`LEO_SCALING_N_SATS` env override)
- `earth_surrogate_gram_scaling.jl` -- GRAM surrogate, same pattern
