# Router evaluation session record — 2026-08-28

Continues `docs/router_evaluation_session_record_20260827.md`. That session found
the Monte Carlo cases had never been verified and fixed the two largest defects.
This one finished the Monte Carlo work, built the instrument that settles it,
and produced the paper's two data tables.

---

## 1. Headline

Monte Carlo regret against the best static route went from **+33..+201%** at the
start of the previous session to **within the measurement floor everywhere**.
Three further defects were found and fixed, and a fourth — the one that would
make the campaign router's wins publishable — was found and deliberately not
fixed, because the benchmark cannot see it.

The adaptive profiles finish at **0.88 of an oracle** that picks the fastest
static route per workload and thread budget after the fact, over the 20
multi-threaded points of the combined 29-point set.

---

## 2. What changed

| commit | change | measured |
|---|---|---|
| `b7ac0a8d` | Hint layer gated on measured work vs measured overhead | mars_aerobraking @8 +18.9% → +2.1% |
| `40a9449b` | Monte Carlo at 1 thread routes to processes | +67..+87% → +1.2..+3.3% |
| `015a18c2` | Hint gate stops re-parsing ENV per decision | predicate 110.2 → 44.5 ns |
| `792cc793` | Campaign bandit: one campaign is one observation, not N | 2.3–6.4× on the shipped campaign path |
| `529ef126` | Outer-split ladder trimmed to plausible widths | router now finds the width optimum |
| `4a361624` | Paired campaign probe | the instrument that settles it |
| `7d387cd7` | Resolution marking on the regret tables | — |

### The campaign router could not learn (`792cc793`)

Two compounding causes. `record_outer_route_feedback!` credited an arm with the
Monte Carlo *sample count*, and the exploration gate read that — so one
64-sample campaign delivered "64 observations" of a route and satisfied any
min-samples guarantee outright. And the first campaign on an arm pays JIT of the
threaded dispatch path and pool spin-up: measured at 0.590 s against serial's
0.060 s, after which every campaign ran serially. One cold sample cost a factor
of 2.6, permanently. `campaigns` now gates exploration while `samples` still
normalises the mean, and the cold timing is evicted once a warm one exists.

### The width ladder offered arms that could only mislead (`529ef126`)

`outer_split_candidates` enumerated `[1, 2, 4, 8, 12]`. `_route_is_proven` stops
once the default has been "tested against **an** alternative", and w1 satisfies
that while telling you nothing. Two arms of five were measured and the ladder
closed at the full budget, where a direct sweep shows width 8 is 47% faster on
`montecarlo_multi_sat` and 55% on `montecarlo_high_accuracy`.

**The obvious fix was wrong.** Forcing the guard to require every candidate be
measured (`SPACEAGORA_OUTER_ROUTE_PROVEN_ALL`, added and default-off) makes all
three campaigns *worse in steady state*, settling on w4 twice. Three noisy arms
can be ranked and five cannot. Trimming the ladder to `[4, 8, 12]` instead lets
the ordinary guard find the optimum on all three.

---

## 3. The instrument

`scripts/paired_campaign_probe.jl` pairs whole campaigns — alternating A,B,A,B
with the order reversed each pair, reduced by a two-sided sign test. It runs the
campaign through the harness's own `ppc_run_sample_batch` rather than
reimplementing it. Validated before use: null control 8/7 at p = 1.0000;
positive control against inner threading 15/0 at −74.5%, p = 0.0001, which
matches the +322% the harness reports for that route.

**Every Monte Carlo residual it was pointed at dissolved.** mars_aerobraking @8
+9.1% → +1.2%; independent_1sat @12 +9.5% → +4.5%; none significant. The one
reading near significance was re-run at 41 pairs and moved away from it.

---

## 4. The resolution floor, and two wrong ways to derive it

**Within-point repeat scatter understates badly and does not predict the
shape.** It is 1.9 points on this data. Robust CV is 4–5× *tighter* at one
thread (0.43–0.63%) than at eight or twelve, which suggests a rung-dependent
band — and direct validation refutes it: at one thread the harness reads +5.8%
and +5.3% where the probe reads +0.5% and −0.7%, errors of 5.3 and 6.0 points,
as large as anything at twelve.

The figure used is the measured disagreement against the paired probes, over
eight comparisons spanning 1, 8 and 12 threads and both workload families:

    0.2, 1.5, 3.6, 5.0, 5.3, 6.0, 7.9, 8.2 points   (median 5.3, max 8.2)

Hence **±8**. The 0.2 is the one comparison where the effect is large (−30.3%
against −30.1%): this harness is accurate where differences are big and
unreliable only where they are small.

**A residual warm-up lands on the first timed repeat** despite the harness
discarding a warm-up solve — at 83% of single-simulation points and **100% of
Monte Carlo points**. It is directional against the adaptive profiles, which pay
one-time costs the statics do not (calibration sweep, hint load): median shift
from dropping repeat 1 is −0.33%/−0.20% for static against −0.70%/−0.59% for
adaptive. Effect on published regret: median −0.14 pp, max 5.31 pp. Medians
absorb it; characterised, not corrected, because discarding a fifth of the data
for a 0.15 pp bias is a bad trade and changing the estimator mid-campaign is
worse. A std-based CV on n=5 reads this warm-up as variance — ~30% on
`multi_128` where MAD reads ~18% and repeats 2–5 read ~11%.

---

## 5. The gap that matters most

**None of the campaign-router work above can appear in any benchmark number.**
`ppc_run_sample_batch` dispatches its own `Threads.@threads` over every sample
and asks the router only for a backend; `select_outer_split!` is called only
from `_campaign_route_plan`, which the harness never enters. Separately,
`ppc_resolve_outer_backend` builds a fresh `OuterRouteState` per point, so the
whole accumulated-statistics half of the router — proven guard, confidence
ranking, width ladder — is unreachable and no change to it can move a published
figure. That property is now documented at the construction site (`3d8afda3`),
along with the warning that persisting state to save startup destroys it
silently.

Closing this gap is the prerequisite for any of the campaign work becoming a
paper result. It is the single highest-value item outstanding.

---

## 6. Open

1. The gap in §5.
2. The harness's residual first-repeat warm-up (§4). Fix deliberately *after*
   the trx50 B9–B14 campaign finishes — not during, or tables get computed two
   ways.
3. `campaign_gnc_6dof_gram` in the paper's capability table is measured at
   `874aea4e`, **39 `src/` files behind** the rest of the table including all of
   `src/parallel/cost/`. It is the only live-GRAM row. A re-run at HEAD matching
   the original configuration exactly is specified and awaiting a decision.
4. `_ppc_apply_cpu_pinning` sizes its mask from thread count alone, ignoring
   `process_workers`, so a process-route point under `--cpu-list` measures the
   pinning artifact. Found by the trx50 session; unfixed.
