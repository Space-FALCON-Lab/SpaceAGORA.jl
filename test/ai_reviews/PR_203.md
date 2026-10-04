# PR201 Follow-up AI Review

## Scope

Independent source review by a separate Codex subagent that did not edit the source patch, against base `eb218881da9c52adeed47256b38fb4231f8afd8c`. Scope is F1 solar ephemeris compatibility/cache selection and F2 malformed thermal CSV handling. Direct heating API policy changes are excluded.

## Changed Files

- `src/simulation/engine/setup.jl`: recognize declared solar requirements while retaining the existing cannonball activation checks.
- `src/vehicle/thermal/thermal_models.jl`: validate each CSV heat value and track grid occupancy independently of numeric values.
- `test/unit/simulation/solar_ephemeris_requirements_tests.jl`: backend, public-entry, initializer/reuse and interpolation regressions.
- `test/unit/dynamics/facet_srp_and_thermal_models_tests.jl`: malformed CSV ordering, value and valid-table controls.
- `test/unit/runtests.jl`: register the solar regression file.

## Findings

No unresolved finding in the reviewed patch. F1 now routes facet and custom declared solar consumers through the existing backend validation and Sun-cache initializer. Cannonball direct/albedo activation remains unchanged; zero-area, IR-only and fully disabled cannonball models retain their exemptions. An inactive first entry does not hide a later solar consumer. Legacy models without a requirement declaration retain the default false capability.

F2 rejects nonfinite/negative heat values before writing and detects duplicates using a separate occupancy mask. Completeness is checked before the matrix is passed to the array constructor, so uninitialized entries cannot reach interpolation. Valid grids and interpolation are unchanged.

## P1 Assessment

No P1 found. The changes affect configuration/cache eligibility and CSV input validation. Force laws, coordinate transformations, thermal formulas and solver equations are unchanged.

Unresolved P1: No

## Tests Added/Updated

The new solar tests exercise actual backend validation and public `run_simulation` rejection, active/inactive/custom/legacy model controls, initialization from the real reuse store, actual midpoint sampling and zero SPICE counters. They restore the touched cache key and environment settings.

CSV regressions cover both NaN duplicate orders, a finite duplicate, standalone NaN/positive infinity/negative infinity/negative values, accepted zero flux, incomplete coverage, and a shuffled valid grid with interpolation. Executed by the primary agent on Julia 1.12.3: 95 focused assertions passed (50 solar and 45 facet/thermal), as did all eight CI shard-planner checks and the public API surface and HPC/extensibility documentation gates. Exact logs are preserved in the accompanying validation evidence.

## Residual Risk

The focused kernel-free tests establish initializer eligibility, cache reuse and sampling; they do not build a fresh cache from native SPICE kernels. The existing native integration suite covers the common cannonball cache builder. No new flight-data, external trajectory or performance campaign is claimed. The pre-existing flat-prefill behavior for inactive cannonball models and optional direct-call NaN API policy are outside these two findings.
