# Router evaluation session record — 2026-08-27 (Monte Carlo)

Follow-on to `docs/router_evaluation_session_record_20260826b.md`. That session
worked the inner axis on single-simulation workloads. This one re-ran the paper
benchmarks, checked the Monte Carlo cases for the first time, and fixed what it
found.

---

## 1. Headline

**The Monte Carlo cases had never been verified, and R4/R5 lost every one of
them.** Regret against the best static route was +33% to +201%. Two defects
account for it, both now fixed, and the remaining regret is 2–23%.

Summed over ten workloads, both adaptive profiles now beat every fixed route:
the best fixed route costs 19% more than R4 and 9% more than R5, and serial
costs 4.5 and 4.1 times as much.

---

## 2. What was wrong

### The Monte Carlo route rule was not a function of per-sample work

`_priority_outer_route_montecarlo` returned `:process` whenever the machine was
medium or large and either the sample count cleared `mc_process_min_samples`
(16) or the mission cleared `mc_process_min_mission_s` (3600 s). Neither term
measures per-sample compute. Sample count rises on both sides of the
threads-vs-process comparison and so cannot discriminate; `mission_time_s` is
*simulated* seconds, and a one-satellite one-hour arc is 3600 by that measure
and 72 ms of actual work. Essentially every campaign of 16+ samples took the
process route.

The process route then lost every time, at 64 samples on 12 threads:

| case | per-sample | process vs threads |
|---|---:|---:|
| `montecarlo_multi_sat` | 0.038 s | 2.45× slower |
| `independent_1sat_1hr` | 0.072 s | 2.90× slower |
| `montecarlo_high_accuracy` | 0.072 s | 1.46× slower |
| `montecarlo_mars_aerobraking` | 0.259 s | 1.41× slower |
| `montecarlo_heavy_aerobraking` | 3.179 s | 1.69× slower |

**There is no crossover in that range.** Even at 3.18 s of compute per sample,
two orders of magnitude above the cheapest case, Distributed's per-campaign
worker startup and per-sample serialisation still cost more than they save. So
no threshold on per-sample cost would have rescued the rule; the premise was
wrong, not the constants.

Native GRAM is unaffected — `_is_native_gram_point_density` decides before this
rule is consulted, and there process isolation is a thread-safety requirement.
`:process` also stays in the *candidate* set so the bandit can rediscover it.

### The telemetry lock was process-global

`PolicyContext` holds a run's telemetry and adaptive state, and
`with_policy_context` gives every `run_simulation` its own. All of it was
protected by one process-global `ReentrantLock`, taken on every policy decision
and every observation — so under a thread-routed campaign all N samples
serialised on one lock, once per RHS call each, guarding state no two shared.

Attribution on `independent_1sat_1hr` against `full_smart`'s +22.4%: the arms
that stop entering the decision path at all both **beat the best static route**
(`nopolicy` −9.4%, recovering 142% of the gap; `innermodes_off` −3.8%, 117%),
which is what identified the bookkeeping rather than any decision it makes.

### The harness was measuring a profile that does not ship

`benchmarks/studies/parallelization_performance/modes.jl` duplicates the routing
knobs `profile_definitions.jl` declares, and the duplicate went stale twice —
R5's scheduler in `622ae2a0` and its `thermal_mode` later. The thermal reversal,
worth −16.2% and −25.3% on two aero workloads, could not reach the benchmark at
all. Fixed, and `ci_ppc_mode_profile_parity_gate` now fails the build on drift
(verified by reintroducing it).

---

## 3. Progression

| case | mc | original | +route fix | +lock fix |
|---|---:|---:|---:|---:|
| `montecarlo_multi_sat` | 16 | +58.3% | −0.6% | **−5.4%** |
| `montecarlo_multi_sat` | 64 | +192.9% | +1.4% | +2.3% |
| `montecarlo_high_accuracy` | 16 | +37.8% | −4.5% | +4.6% |
| `montecarlo_high_accuracy` | 64 | +43.4% | −0.8% | **−0.7%** |
| `montecarlo_mars_aerobraking` | 16 | +45.2% | +9.2% | +9.4% |
| `montecarlo_mars_aerobraking` | 64 | +41.3% | +13.2% | +3.3% |
| `independent_1sat_1hr` | 16 | +171.1% | +17.6% | +8.9% |
| `independent_1sat_1hr` | 64 | +201.0% | +22.8% | +18.7% |

---

## 4. Consolidated result, ten workloads, twelve threads

Best static varies: R1_a on 3, R3 on 4, R2 on 2, R1_b on 1. Worst fixed route
runs 1.3× to 6.2× the best on a given workload.

| profile | total (s) | overhead |
|---|---:|---:|
| R4 | 22.59 | +0% |
| R5 | 24.74 | +10% |
| R3 | 26.90 | +19% |
| R1_a | 27.58 | +22% |
| R1_b | 50.48 | +123% |
| R2 | 75.01 | +232% |
| R0 | 102.04 | +352% |

Per workload R5 is faster than the best fixed route on 3, within 2% on 2, and
slower on 5 (by 2–23%). **It is not a uniform per-workload win**, and the
section in the paper says so. The aggregate win comes from the three workloads
where inner routing is decisive (30–67% faster than any fixed route) against a
bounded loss elsewhere.

All modes pass trajectory parity, `full_smart` included.

---

## 5. Open: R4 beats R5, and the hint layer is why

R4 and R5 are within 1% of each other on the three decisive wins, so R5's
persistent hints and measured-reward selection contribute nothing where the
gains are. Disabling the hint layer and changing nothing else:

- `montecarlo_heavy_aerobraking` −12.0% (essentially the whole R5−R4 gap)
- `gravity_4096sat_l50_vacuum_1hr` −5.0%
- `independent_1sat_1hr` **+12.8%** — the layer helps here

Two cases better without, one clearly worse. Not resolved, and deliberately not
acted on: choosing a side on 2-against-1 would be exactly the unmeasured profile
constant this codebase has now had to revert twice (R5's scheduler, R5's
thermal mode). **R4 is the better default on the workloads measured.**

---

## 6. Method notes

- The harness runs modes in blocks with a fresh process per point, which the
  2026-08-26 record warns fakes effects. It is appropriate here only because the
  Monte Carlo differences are 1.4× to 3×, far outside that noise. Anything
  marginal still needs `scripts/paired_profile_probe.jl`.
- `make_router_regret_table.py` keyed on `(case, threads, mode)` and so took a
  median *across* the Monte Carlo sample ladder, reporting one number for a
  16-sample and a 64-sample campaign. `mc_samples` is now part of the key.
- GRAMSuite is a *weak* dependency: it is in `[weakdeps]`/`[extras]` with a
  `[sources]` path and not in `[deps]`, so `using GRAMSuite --project=.` fails
  by construction. That is not a broken install. The harness loads it by pushing
  the vendored path onto `LOAD_PATH` (`cases.jl:23`), which attaches
  `SpaceAGORAGRAMSuiteExt`. The native `libGRAM.so` is built and the submodule
  env is instantiated.

## 7. Not done

1. GRAM-backed Monte Carlo (`montecarlo_mars_gram_live`, `campaign_gnc_6dof_gram`)
   was not re-run. Those are the cases where the process route is a requirement,
   so the reversal above should not touch them, but that is reasoning rather
   than measurement.
2. The thread ladder was fixed at twelve. The earlier consolidated table swept
   one, eight and twelve; this one does not.
3. `montecarlo_heavy_aerobraking` and `independent_1sat_1hr` remain ~20% behind
   the best fixed route under R5 (~0–1% under R4). See §5.
