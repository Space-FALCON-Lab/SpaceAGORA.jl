# Router evaluation session record — 2026-08-25/26

Working record for the adaptive-policy work on `router-eval-expanded-b6`. Picks
up from `docs/router_evaluation_session_record_20260824.md`.

Read the methodology section before trusting any number in here or generating
new ones. Several results in this session were reported and then retracted, and
in every case the cause was the measurement rather than the code.

---

## 1. Headline

R4/R5 now carry **no statistically significant regression against the best
static route** on any workload measured with an instrument able to resolve one.

Final state, `full_smart` against `inner_only`, paired probe, 21 pairs, 12 threads:

| case | median A/B | pair wins | p | verdict |
|---|---:|---|---:|---|
| `gravity_4096_l50` | +5.4% | 8 / 13 | 0.383 | not distinguishable |
| `heavy_1024_l50_6hr` | +0.4% | 8 / 13 | 0.383 | not distinguishable |
| `interact_256_aero` | **−31.3%** | 21 / 0 | 0.0000 | **R5 faster** |

`gravity_4096_l50` was +7.8% at p = 0.0002 (2 wins against 19) before the final
fix. It is now indistinguishable from zero.

---

## 2. What changed, in dependency order

| commit | change |
|---|---|
| `3b7aac25` | Forced-region short-circuit extended; inner scheduler made a routed choice |
| `ba8346e8` | Flat-queue partials reduction parallelised — up to **5.3× at wide splits**, bit-identical |
| `9d615cc7` | Outer-route bandit persisted across processes; proven-default guard |
| `209bc3dc` | mc=1 inner-policy regression: **+84% → −6.5%** |
| `d4f83343` | Successive-halving sweep (200 → 99 RHS calls); multibody fast path |
| `ed81ad66` | No-regret floor on plan selection; calibration signature v2 → v3 |
| `e6dc7539` | Outer/inner thread split adapted rather than computed |
| `f1bc92aa` | `parallel_cost.jl` allowlisted; `pr_runtests.jl` green end-to-end |
| `f74e4199` | Paired in-process probe |
| `622ae2a0` | R5 defaults to the static scheduler — last significant regression eliminated |

Cost-model commits (`2d939eef` … `0b457774`, `ef9e298f`, `c06f2ed4`) are
described in §5; that work is **not shipped**.

---

## 3. Measurement methodology — read this first

Four independent failures this session, each of which produced a wrong number
that looked plausible.

**Block ordering fakes effects.** A sequential A/B on the thermal callback gave
**+11.1%**; the same comparison interleaved and paired gave **+1.7%**. The 11.1%
was entirely measurement order. Never compare profiles in blocks. Use
`scripts/paired_profile_probe.jl`, which alternates A/B/A/B *and* alternates the
order within each pair.

**Per-solve spread is 15–45%; the residuals under investigation are 3–10%.** More
repeats does not fix this — the variance is per-solve, not per-session. Pairing
does fix it, because drift is common-mode within a pair.

**R5 needs more warm-up than the static profiles.** It pays a calibration sweep
and hint-state load they do not. Cutting `--warmup` from 2 to 1 let its setup
leak into its timed samples and inflated its apparent regret by ~3 points. Any
comparison involving R5 needs `warmup ≥ 2`.

**A cached calibration entry silently bypasses a fix to calibration.** A cache
hit returns before the sweep runs, so the no-regret floor never executed and the
first re-measurement of it tested nothing. Any change to sweep *behaviour* needs
the signature version bumped, or it will not reach an already-calibrated machine.

Related: report **regret**, not decision accuracy, wherever there is a choice.
Accuracy is binary and undefined-ish at near-ties, so it swings run to run;
regret is continuous and degrades gracefully.

---

## 4. Tooling

- `scripts/paired_profile_probe.jl` — the instrument that settles marginal
  differences. One process, interleaved, sign test, explicit "not
  distinguishable" verdict. Refuses fewer than 6 pairs (a sign test on n pairs
  bottoms out at 2/2ⁿ, so 5 can never clear 0.05). Derives its env from
  `ParallelProfiles.profile_env_pairs`, which is the shipped source of truth —
  an earlier duplicated table drifted within an hour of being written.
- `scripts/make_router_regret_table.py` — consolidated LaTeX regret table from
  paper-benchmark raw CSVs. Regret is against the best static route *observed*,
  not a nominated baseline.
- `scripts/calibrate_machine.jl`, `scripts/validate_cost_model.jl` — cost-model
  calibration and cold validation. See §5.

The benchmark harness spawns a **fresh Julia process per (case, mode, thread)
point** — ~40 s startup against a few seconds of solving. For anything marginal,
prefer the paired probe. A useful lighter harness configuration is 3 modes ×
9 repeats rather than 6 × 3: half the wall time, three times the samples.

Known harness defect (pre-existing): `--parity-cases=` does not restrict the
parity list; an empty value falls back to the default cases.

---

## 5. Cost model — built, validated, NOT shipped

`src/parallel/cost/`, driven by `scripts/validate_cost_model.jl`. `select_plan`
is called by no routing path.

Decision accuracy **30–45%** across six runs on identical constants
(30/30/35/35/35/45) — quote it as a range. Predictions are deterministic given
cached constants, so the spread is measurement noise reordering near-ties.

It is useless as a chooser and usable as a **filter**: the true best sits in its
top-2 65% of the time and top-3 85%, both at **0.0% median regret**. Worst rank
of truth was 8 of 10 candidates, so a hard cap on a swept set must be chosen
against that tail, not the median.

Six modelling errors found by cold validation, each of which looked plausible
until measured: scope mismatch (whole-RHS vs effector pass), per-mechanism
dispatch (Polyester 704 ns vs channel pool 18.4 µs at 12 workers), flops-per-lane
vs loop passes (21×), a latency-bound calibration kernel (24×), lane rate indexed
by workspace footprint not batch width, and assumed-linear parallel scaling
(real speedup saturates near 5×).

**Hold the real kernel out of the fit.** `simd_terms` and `coeff_touches` stay
0.998 correlated across every configuration the real kernel can produce (exactly
1.000 zonal), so a joint fit on real workloads cannot separate them — and two
splits predicting `flat` identically can predict `satellite_batch` 49% apart.

---

## 6. Claims made and later retracted

Recorded because the reasoning that produced them was plausible.

- *"The N×W partials reduction is what breaks the L50 cases."* Wrong. A
  harmonics-only constellation returns at `_count_flat_queue_only_effectors == 0`
  before the reduction runs. The reduction fix is real and valuable, but for
  multi-effector queues.
- *"The vacuum-gravity losses are consultation overhead."* Wrong. Telemetry shows
  `policy_decisions_total = 0` for *every* profile on those workloads.
- *"`thermal_mode=on` costs the gap."* Artifact of block ordering (11.1% → 1.7%).
- *"The no-regret floor made things worse."* The measurement was invalid — a
  cached entry bypassed the sweep entirely.
- *"Decision accuracy is 30–35%."* Understated; it is 30–45%.

---

## 7. Open items

1. **Paper `\input` is still commented.** `../JAIS-2026-SpaceAGORA/main.tex`
   §VII.A carries the consolidated subsection and a commented
   `\input{data/router_regret_consolidated.tex}`. The generated table is built
   from the **pre-fix** 33-row run and no longer describes the shipped router.
   Either re-run the broad set (~7.5 h) or ship the paired-probe numbers for the
   three cases with a note that the full sweep predates `622ae2a0`.
   §IV was updated in place and does describe the current implementation.
2. **Route persistence and split adaptation are unmeasured.**
   `ppc_resolve_outer_backend` constructs a fresh `OuterRouteState()` per
   measurement point *by design*, so neither can appear in any harness number.
   They act only on repeat runs with accumulated state. Unit tests cover the
   mechanisms; wall-clock value is unproven.
3. **Pre-existing test failure**, unrelated to this work and verified against a
   stashed baseline: "Adaptive Campaign Routing" in
   `test/suites/02_callbacks_parallel_and_smoke_tests.jl` errors with
   `UndefVarError: SimulationCampaigns / ParallelProcess not defined in Main` on a
   Distributed worker. Suites are raw-included into the coordinator's `Main`, so
   closures sent to workers reference names the workers never receive.
4. **`ci_canonical_path_contract_gate` trips on untracked local files.** It walks
   `SCAN_ROOTS`, so any scratch doc under `docs/` with forbidden path tokens
   fails it — currently `docs/IMPLEMENTATION_SUMMARY.md`. Environmental, not a
   code defect. Move the file aside to run `pr_runtests.jl` clean.
5. **`gravity_4096_l50` still shows +5.4% median** at 8/13 pairs. Not
   significant, but not proven zero either. If it matters, more pairs.

---

## 8. Coordination

A peer session (`spaceagora-jl-56`) is working on an inner-policy cost hierarchy
plan at `docs/architecture/inner_policy_cost_hierarchy_plan.md`. It found the
`parallel_cost.jl` allowlist gap. Agreed boundaries: it does not edit
`src/parallel/policy/adaptive_decision.jl` or `src/parallel/cost/`.

Its plan's L0 gate must sit at or above the forced-region short-circuits, be
cheaper than ~9 ENV reads plus two lock acquisitions, and stay
measurement-transparent — a short-circuit may skip the decision, it may **not**
skip the measurement feeding the next decision. `heavy_work` comes from
`_update_effector_cost_model!`, gated by `effector_decision.policy_applied`, so a
fast path that skips the measurement starves the signal it gates on and latches
silently and permanently.

---

## 9. Machine hygiene

Killing a benchmark coordinator does **not** kill its Julia workers — they
reparent and keep burning a core, which then contaminates whatever runs next.
This bit once in this session. Always verify:

```
pgrep -af julia | grep -v languageserver
ps -eo pcpu,etime,comm --sort=-pcpu | head
```
