# Inclination Experiment: 50 Helpers, 50 Orbits, One-Degree Target

For repeated first-call and warmed-up wall-clock measurements at 50 and
100 orbits, see the [runtime comparison](RUNTIME_COMPARISON_50_100_ORBITS.md).

## Setup and implementation availability

This repeats the earlier 50-helper experiment with the **inclination-increasing
scheduler `:gve_inc`** and initial **target inclination of 1 degree**.

- Target: initial altitude **1050 km**, eccentricity 0, inclination 1 degree.
- Helpers: 50 initially circular, equatorial orbits at **1000 km**.
- Helper phases: 0, 7.2, ..., 352.8 degrees.
- Target RAAN, argument of periapsis, and initial true anomaly: 0 degrees.
  The target therefore begins at its ascending node.
- Duration: 50 initial target orbital periods,
  **318,567.018024937 simulation seconds**.
- Other settings unchanged: 227 kg each, 10 kW, magnification 100,
  200 km range, central gravity, no drag or J2.
- Earth radius: 6,378,136.6 m; propagation gravitational parameter:
  3.98600436233e14 m^3/s^2.
- v1 uses `eta=1`, `beta=1` so force magnitudes match the prototype.

| Requested implementation | Outcome |
|---|---|
| Prototype | Ran its existing `:gve_inc` scheduler; `Success` |
| 4 / v1 | Ran its existing `:gve_inc` scheduler; `Success` |
| 9 integration | Ran the newly added physical-rate `:gve_inc` policy; `Success` |

The prototype/v1 results are retained from the earlier experiment. The
current workspace adds 9's inclination policy and its new run to this
comparison; no SMA-policy substitution was used. See the
[source comparison](INCLINATION_SCHEDULER_COMPARISON.md) for formulas and
the explicitly defined equatorial boundary behavior.

## Results

| Metric | Prototype | v1 | 9 |
|---|---:|---:|---:|
| Initial target inclination (deg) | 1.0000000000 | 1.0000000000 | 1.0000000000 |
| Final target inclination (deg) | **1.0021997427** | **1.0021991073** | **1.0022020161** |
| Inclination increase (deg) | **0.0021997427** | **0.0021991073** | **0.0022020161** |
| Inclination increase (arcseconds) | 7.919074 | 7.916786 | 7.927258 |
| Target SMA change (m) | +0.524977 | -1.258155 | -2.703304 |
| Final eccentricity | 0.000006938538 | 0.000006991004 | 0.000007055112 |
| Recorded active segments | 50 | 50 | 50 |
| Startup transients under the analysis definition | 0 | 0 | 0 |
| Link active at simulation end? | No | No | No |

The inclination gains differ by **0.000000635453 degrees**, approximately
**0.002288 arcseconds**, or **0.028892% of their mean gain**.

Across all three, the largest inclination-gain spread is
**0.000002908861 degrees**, or **0.132204% of the mean gain**.

All three schedulers therefore produced very similar inclination increases in
this scenario. Their final SMA changes differ in sign, but this policy
optimizes inclination, not SMA. The approximately 1.783 m SMA difference
between prototype and v1 (3.228 m across all three) is not evidence that a
scheduler is trying to lower SMA; SMA is a
secondary outcome of the forces and numerical integration.

## Helper selections and timing

After converting raw spacecraft indices to common physical helper slots,
all three runs have the same 50-segment helper sequence:

```text
1,
50, 50, 49, 49, 48, 48, 47, 47, 46, 46, 45, 45, 44, 44,
43, 43, 42, 42, 41, 41, 40, 40, 39,
38, 38, 37, 37, 36, 36, 35, 35, 34, 34, 33, 33, 32, 32,
31, 31, 30, 30, 29, 29, 28, 28, 27, 27, 26, 26
```

Repeated helper numbers mean separate active intervals with a gap, not
simultaneous duplicate links. The phase of physical slot $j$ is
$7.2(j-1)$ degrees. The prototype uses raw helper index $j$, while v1 uses
raw index $j+1$, as does 9, because their target occupies spacecraft index 1.

| Timing diagnostic | Prototype | v1 |
|---|---:|---:|
| First contact: slot 1 starts (s) | 0.0001420725 | 0.0000124336 |
| First contact ends (s) | 1,592.923232 | 1,602.311103 |
| Last contact: slot 26 starts (s) | 315,358.884413 | 315,362.825261 |
| Last contact ends (s) | 315,457.666326 | 315,464.433449 |
| Estimated active time (s) | 60,067.363734 | 60,000.007531 |
| Estimated duty fraction | 18.8555% | 18.8343% |

9's first slot-1 contact starts at 0.0000111074 s and ends at
1,595.093641 s. Its final slot-26 contact spans
315,362.186373 to 315,472.186373 s. Recorded active time is
60,102.186362 s, giving an estimated duty fraction of 18.8664%.

All these timestamps and durations are **simulation time**.
The first contact starts very early but persists for about 1,600 seconds:
it is a real contact, **not** a microsecond-long startup transient.
The runner's earlier diagnostic classified activations by start time alone;
that was corrected to distinguish a brief initialization selection from a
genuine early contact.

For prototype versus v1, the largest observed activation-time difference
is **13.410573 s**, and the largest observed deactivation-time difference
is **9.958556 s**.

Holding each recorded selection until the next timestamp gives
**99.847017% time-weighted schedule agreement** over the whole run,
including inactive intervals. The corresponding disagreement time is
487.353303 s. This is agreement between recorded schedules, not an exact
comparison of every internal RHS force evaluation.

The additional comparisons, with the same measurement method, are:

| Pair | Recorded-schedule agreement | Maximum start difference (s) | Maximum end difference (s) |
|---|---:|---:|---:|
| Prototype vs 9 | 99.810196% | 21.515271 | 17.099718 |
| v1 vs 9 | 99.888170% | 11.223785 | 11.863506 |

## Interpretation and caveats

The prototype and v1 rank directional sensitivity:

$$
q_i=\frac{r\cos u}{h}
\left(\hat{\mathbf{f}}\cdot\hat{\mathbf{n}}\right),
\qquad
\dot{i}=\frac{F_{\mathrm{mag}}}{m}q_i.
$$

9 instead ranks the actual physical rate $\dot{i}$ computed from angular
momentum and its derivative, including force/mass scaling. Here equal link
strengths and one fixed target make the positive scaling common to all
candidates at a given state. The observed agreement is therefore expected;
it does not prove identical behavior for unequal powers or multiple targets.

All select strictly positive scores. Since the helpers remain
near the equatorial plane while the target is inclined, link directions
can have a component normal to the target's orbital plane. This provides
inclination-changing opportunities that the previous coplanar setup lacked.
The sign and size of the contribution also depend on argument of latitude;
being in range does not mean the link should stay active for the whole
encounter.

- Prototype: its original Vern9 path uses `reltol=1e-12`, `abstol=1e-12`,
  and performs selection during RHS force evaluation. Its schedule CSV is
  reconstructed at 6,349 saved states, not at every RHS evaluation.
- v1: the existing wrapper executes actual v1 source with a central-gravity
  ODE, Tsit5, `reltol=1e-9`, `abstol=1e-11`, and a 10 s maximum step.
  It records 32,028 accepted-step scheduler updates.
- 9: native `run_simulation` uses Tsit5, the same requested orbit tolerances
  and 10 s step cap, and records 31,868 scheduler updates. Its median and
  maximum update intervals are 10 s. Measured simulation-call runtime was
  11.271 wall-clock seconds, not including package loading or result analysis.
- Prototype schedule spacing has a median of 20.226 s and maximum of
  124.448 s; v1's median and maximum are 10 s.
- Availability rules and the v1 scorer's hard-coded gravitational parameter
  retain the differences described in the earlier comparison.
- 9 uses its existing maximum-total-score terminal matching. Only
  target-helper pairs are registered and the target has one terminal, so it
  chooses at most one helper, like the other two in this experiment.
- No timestep-convergence study was performed. Similar outputs do not
  establish exact trajectories or a numerical error bound.

## Reproduce and inspect

Run from the workspace root using the existing Julia projects:

```sh
GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/run_1_prototype.jl 50 50 gve_inc 1
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/run_4_v1.jl 50 50 gve_inc 1
GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project="9_SpaceAGORA.jl-main_full integration" SciTech_Paper/run_9_full.jl 50 50 gve_inc 1
julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/compare_sma_runs.jl 50 50 gve_inc 1
```

Julia version: 1.12.1. Despite its original SMA-oriented name, the
[comparison script](compare_sma_runs.jl) now also supports this three-code
inclination case. It validates the requested policy, inclination, altitudes,
duration, return codes, summary arithmetic, schedule shape, helper mapping,
and contact counts before writing derived statistics.

Results are separate from all prior SMA runs:

- Prototype: [summary](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/prototype/summary.toml),
  [schedule](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/prototype/schedule.csv).
- v1: [summary](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/4_v1/summary.toml),
  [schedule](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/4_v1/schedule.csv).
- 9: [summary](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/9_full/summary.toml),
  [schedule](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/9_full/schedule.csv).
- Derived [comparison statistics](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/comparison.toml)
  and [all contact intervals](results/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/contacts.csv).

The new runner arguments default to the original SMA policy and zero target
inclination when omitted; previous commands remain supported.

### Implementation validation

The new standalone
[inclination policy tests](../9_SpaceAGORA.jl-main_full%20integration/test/unit/simulation/inclination_scheduler_tests.jl)
pass 232 checks covering GVE equivalence, finite differences, force/mass
scaling, equatorial one-sided rates, invalid inputs, endpoint selection,
availability, positive-score matching, and shared SMA-scoring regression.

```sh
julia --project="9_SpaceAGORA.jl-main_full integration" "9_SpaceAGORA.jl-main_full integration/test/unit/simulation/inclination_scheduler_tests.jl"
```

The existing broad interlink test file has a pre-existing include of a
missing `test/examples/gve_sma_interlinks.jl`; it was not used as a passing
validation claim. The standalone tests and the full 50-orbit native run
validate the new policy independently of that broken test prerequisite.
