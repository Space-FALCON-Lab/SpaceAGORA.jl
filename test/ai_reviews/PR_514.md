# Review record: process-local cache reset

## Scope

Repair a precompile-state snapshot failure reproduced on unchanged main
`0f03bbfcef0aa694fc01cd3209fed5aa07edafd7`. Author verification, not independent review.

## Changed Files

- `src/precompile_workload.jl`: reset planet-frame ephemeris reuse caches through
  their locked helper and clear the multibody-mode cache.
- `src/orchestration/simulation/campaigns/predictive_planner.jl`: reset predictive
  constants generation with its values under the existing lock.
- `test/unit/parallel/precompile_state_tests.jl`: exclude the native-lock entry
  site scratch from snapshot equality alongside acquisition time and wait time.

## Findings

The unchanged-main snapshot reported five differences after workload/reset:
native-lock acquisition site, predictive-constants generation, ephemeris reuse
cache, and two names for the same multibody-mode cache. The reset now returns the
three caches to their unread state. The lock site is overwritten at each tracked
acquisition before use, so it is transient scratch rather than a cached host fact.
The actual caches remain in the snapshot; the repair clears them in production.

## P1 Assessment

No P1 findings identified in the changed reset paths. No numerical model changed.

## Tests Added/Updated

Julia 1.12.1: unchanged main fails the first snapshot testset (4 pass, 1 fail).
The repair passes all 5 snapshot assertions and all 8 direct reset assertions.
All 33 PR contract gates pass. The full unit tree is still running; hosted CI,
coverage and release checks are pending. No full-suite pass is claimed.

## Residual Risk

The reset is used during precompilation and test setup. Full integration checks
must finish before merge. The older checkpoint/output branches exposed this
baseline failure; their full unit runs must not be called passing because this
separate repair passes its focused checks.

Unresolved P1: No
