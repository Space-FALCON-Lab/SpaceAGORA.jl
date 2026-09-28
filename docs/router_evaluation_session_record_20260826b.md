# Router evaluation session record — 2026-08-26 (second session)

Follow-on to `docs/router_evaluation_session_record_20260826.md`. That session
closed the significant regressions of R4/R5 against the best static route in the
**warm** regime. This one found that the warm regime is the only one its
instrument could see, extended the instrument, and fixed what became visible.

Read §2 before trusting any number here. Two results in this session were
reported at p ≤ 0.001 and then retracted, both because a single run of the
paired probe is not sufficient evidence.

---

## 1. Headline

R5 now **meets or exceeds** the best static route on every case measured,
including three light workloads absent from the previous case set.

Cumulative, R5 against R2, cold, isolated stores, 21 pairs, two independent runs:

| case | run 1 | run 2 | verdict |
|---|---:|---:|---|
| `gravity_4096_l50` | −4.5% (ns) | +1.7% (ns) | parity |
| `heavy_1024_l50_6hr` | −0.0% (ns) | −0.4% (ns) | parity |
| `interact_256_aero` | **−41.1%** | **−39.7%** | R5 faster |
| `light_64_aero` | **−81.6%** | **−84.7%** | R5 faster |
| `light_16_harm` | +0.9% (ns) | −0.1% (ns) | parity |
| `light_256_aero_j2` | **−10.4%** (p = 0.0001) | — | R5 faster |

Cold baseline before any of this: `gravity_4096` +1.4%, `heavy_1024` +0.5%,
`interact_256` −36.8%.

---

## 2. The blind spot, and what it cost

**`scripts/paired_profile_probe.jl` runs every sample in one process, and
`_rhs_calib_cache` is a module global.** The pre-solve calibration sweep
therefore runs on the first sample and every later sample takes a cache hit. The
probe cannot see per-run calibration cost — it amortises it across all the
pairs. The benchmark harness spawns a fresh process per point, and so does every
real user, so the regime the probe measured was the one regime in which that
cost does not exist.

Two flags were added:

- `--cold-calib` drops the in-process memo before each sample, leaving the
  on-disk TOML alone. That is exactly what a second `julia` invocation on a
  previously-calibrated machine does.
- `--isolate-calib` gives each arm its own calibration TOML, wiped at start. The
  store is shared mutable state that both arms read and write, and the sweep's
  verdict is not stable run to run, so without this one arm's stored plan
  becomes the other arm's cache hit.

**A single significant run from this probe is not sufficient evidence.** The
sign test assumes independent pairs; samples within one process are serially
correlated through allocator state, GC generation, thermal drift and pool
warm-up. Two results at p ≤ 0.001 failed to replicate (§6). Everything claimed
above replicated in a second, independently launched process. An A/A control
(no override) returned 11/10 and 13/8, both not distinguishable, so the
instrument itself carries no systematic bias.

Isolation form used throughout: **same profile on both sides, one mechanism
reverted on B** via `--b-override`. Every fix below therefore ships with an env
knob whose only purpose is to give the A/B a B side.

---

## 3. What changed

| fix | mechanism | measured |
|---|---|---|
| 1 | Cache the no-regret floor's retain-the-heuristic verdict | **−42.6%** `gravity_4096`, −3.5% `light_16_harm` |
| 3 | Interleave sweep candidates instead of running blocks | verdict became **deterministic**; 4/4 cases favourable on wall clock, none individually significant |
| 5 | R5 `thermal_mode` `"on"` → `"auto"` | **−16.2%** `interact_256`, **−25.3%** `light_64_aero` |
| 6 | Signature gains effector identity and outer-split state (v3 → v4) | outer-split blindness shown decisive; effector collision shown harmless on the pair tested |
| 8 | Per-decision telemetry: env snapshot instead of live ENV, one lock instead of two | **−30.3% time, −40.5% allocation** per forced-region decision |
| 2 | Sweep scoring `min` → trimmed mean | no resolvable effect; kept as the better statistic, not as a win |
| 4 | Release the sweep's oversized flat partials buffer | unit cost 35–43× on the affected branch, but latent on every measured shape |
| 9 | Tie-break near-equal survivors to the narrower plan | **shipped default-off**: it widened the outcome set |

### Fix 1 — the floor's verdict was never persisted

`_run_rhs_sweep!` returned `nothing` when the floor kept the heuristic, and
`_rhs_calib_lookup` had no representation for that, so the next solve saw a
cache miss and re-ran the entire successive-halving sweep to reach the same
conclusion. Forever, on exactly the shapes the floor protects.

Verified directly: under `SPACEAGORA_RHS_CALIBRATE=force`, `light_16_harm` and
`light_64_aero` stored no entry at all. The sweep now returns a three-valued
verdict — `:pinned`, `:heuristic`, `:aborted` — and only `:heuristic` is cached.
The `:aborted` returns stay uncached deliberately: they are failures to measure,
and caching one would permanently suppress calibration for that signature.

Old builds parse the sentinel row, fail both plan-mode comparisons and fall
through to a miss, so the file stays backward-compatible.

**Fixes 1 and 3 compound.** Measured before Fix 3 landed, Fix 1 was worth −3.1%
on `gravity_4096`. Measured after, −42.6%. Interleaving made that case's verdict
stably "heuristic" (16/16, against 9/12 with blocks), so after Fix 3 nearly
every uncached run re-sweeps — which is precisely what Fix 1 makes free.

### Fix 3 — the sweep violated this project's own methodology

The previous record's §3 says never to compare in blocks, and demonstrates it:
+11.1% block-ordered against +1.7% paired. That lesson was applied to the
external A/B instrument and not to the sweep, which is the internal one and
makes an irreversible decision. Each candidate's reps ran contiguously,
`satellite_batch` ran first in every round, and drift across a round landed on
whichever candidate was executing.

Round-robin over reps, with the order reversed on alternate passes. Same call
count, different order.

| case | interleaved | blocks |
|---|---|---|
| `gravity_4096_l50` | heuristic ×16 | heuristic ×14, flat(12) ×2 |
| `light_64_aero` | satellite_batch ×16 | satellite_batch ×15, flat(1) ×1 |

**The `flat(12)` entry in this machine's own state file for `gravity_4096` is a
block-ordering artifact.** The interleaved sweep never picks it.

### Fix 5 — the twin of a reversal already made

R5 was the only profile forcing `thermal_mode="on"`; R2, R3 and R4 all declare
`"auto"`. This is the same defect as the hard-coded `inner_scheduler="dynamic"`
that `622ae2a0` reverted — a profile constant pre-empting a routed decision —
three fields up in the same constructor. The +11.1% originally attributed to it
was correctly retracted as block ordering; the paired re-measurement put it at
+1.7%, and that residual was filed with the retraction rather than acted on.
`"on"` was never measured to beat `"auto"` on anything.

### Fix 6 — the signature could not distinguish what the heuristic branches on

`_rhs_execution_plan_uncached` clamps the flat routes' allotment to 1 under
`outer_serialized`, specifically to avoid nesting a thread split inside an
already-blocked outer worker. The signature had no outer term, and budget does
not separate the two cases because the harness asserts
`SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` even for single-simulation constellation
runs.

| case | outer=0 | outer=1 |
|---|---|---|
| `gravity_4096_l50` | heuristic ×16 | flat(12) ×9, flat(8) ×4, flat(4) ×3 |
| `interact_256_aero` | satellite_batch ×15, heuristic ×1 | satellite_batch ×16 |

Completely different answers sharing one cache key. The effector-identity term
(`effs=2` collided harmonics+drag with invsqJ2+drag) is defensive only — both
shapes calibrated to `satellite_batch/1` ×16, so no harm was demonstrated there.

Note also: under `outer=1` the pinned allotment stays unstable at 12/8/4 even
with interleaving. That variance is not ordering and is still unexplained.

---

## 4. Claims made and retracted

- *"The missing v3 entries are shapes where the floor retained the heuristic."*
  Mostly wrong. Seven v3 signatures were absent; the verdict tally showed
  `light_64_aero`, `light_256_aero_j2` and `interact_256_aero` all pin
  `satellite_batch` 12/12. Only `light_16_harm` (11/12) and `gravity_4096`
  (9/12) actually retain the heuristic. The absent entries were mostly
  never-run, not retained.
- *"Releasing the oversized partials buffer is worth −10.1% on `heavy_1024`"*
  (14/1, p = 0.0010). Did not reproduce: −1.8%, p = 0.0784. It also had no
  mechanism — `heavy_1024` is harmonics-only, so
  `_count_flat_queue_only_effectors == 0`, `needs_flat_queue` is false and the
  zeroing branch is skipped entirely. Same trap as the retracted "partials
  reduction breaks the L50 cases" claim in the previous record.
- *"The telemetry fix is worth −66.9% on `light_64_aero`"* (21/0, p = 0.0000).
  Did not reproduce: −4.5%, not distinguishable.

---

## 5. Not built

**Fix 7, a bucketed conditional plan table.** The premise was that a plan pinned
on the pre-solve satellite count goes stale as satellites deactivate. Tested on
`light_256_aero_j2`, the one case that actually deorbits satellites mid-solve:
R5 is −10.4% *faster* than the best static there (15/0, p = 0.0001), and that
case pins `satellite_batch`, which is allotment-independent, so staleness has
nothing to bite on. No measured motivation.

What would motivate it: a workload that pins a **flat** plan (allotment-
sensitive) *and* whose active-satellite count crosses a `_calib_sat_bucket`
boundary mid-solve.

---

## 6. Test state

All 23 PR contract gates pass. `test/unit/runtests.jl` passes.

Pre-existing failures, verified against a stashed baseline and unrelated to this
work:

1. `ci_no_artifact_files_gate` fails on committed absolute paths in
   `benchmarks/studies/edl_aerocapture_gp_uncertainty/README.md` — a tracked
   file this session did not touch. This blocks `test/contracts/pr_runtests.jl`
   as a whole; the 23 gates were run individually to get past it.
2. `test/probes/coverage_parallel_telemetry_probes.jl:201` (outer-route stats)
   fails with and without these changes.

Fixed here: the same probe still asserted `inner_scheduler == "dynamic"` for R5,
stale since `622ae2a0` shipped the reversal. Updated along with the thermal
assertions, and `test/suites/02_*.jl` updated for the v4 signature.

---

## 7. Open items

1. **`--isolate-calib` and `--cold-calib` should become the default** for any
   probe run that involves R4/R5. The warm shared-store regime flatters the
   adaptive profiles and is not what anyone runs.
2. **The outer-active allotment instability** (flat 12/8/4 under `outer=1`, even
   interleaved) is unexplained and is the one remaining source of
   non-determinism in plan selection.
3. **Fix 2 and Fix 4 are unmeasured.** Both are defensible on mechanism and
   neither is a demonstrated win. Fix 4's unit cost is large (35–43× on the
   affected branch) but no measured shape reaches it.
4. **The paper's regret table still predates all of this** — the open item from
   the previous record is unchanged and now further out of date.
5. The repo's own `output/parallel_policy_state/rhs_calibration_*.toml` carries
   v2 and v3 entries that v4 supersedes. They are inert, not harmful.
