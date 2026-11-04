# Five-Orbit SMA Scheduler Comparison

## Setup

This comparison uses 20 circular, equatorial helpers at 1000 km altitude and
one circular, equatorial target at 1050 km altitude. The target starts at true
anomaly 0 degrees; helpers start at 0, 18, ..., 342 degrees. The run ends after
five target orbital periods.

All runs use 227 kg spacecraft, a 200 km link limit, 10 kW input power,
magnification 100, central point-mass gravity, and no drag or J2. The shared
propagation constants are Earth radius 6,378,136.6 m and gravitational
parameter 3.98600436233e14 m^3/s^2. Laser efficiency and reflectivity are 1 in
the v1 model, matching the prototype and 9 model force magnitudes.

The 9 and v1 integrations use Tsit5 with a 10 s maximum step. The prototype
uses the Vern9 integration path in its existing runner; its saved solution has
fewer points, while its SMA scheduler is evaluated from within the RHS.

## Results

All three integrations returned `Success`. All selected the same physical
spacecraft: helper slot 20, initially at true anomaly 342 degrees. The link
was active only near the end of the run, at roughly 4.502 to 4.912 target
orbits.

| Implementation | Final target SMA change | Scheduler/link observations | Link activity |
|---|---:|---|---|
| Kuang prototype | +59.89297 m | 350 saved states; 1 activation and 1 deactivation observed by rescoring those states | Sampled interval 28,681.10–31,284.50 s |
| 4 / v1 | +59.66245 m | 3,196 accepted-step scheduler updates; 1 physical activation | 28,689.10–31,289.10 s |
| 9 integration | +59.83371 m | 3,196 accepted-step scheduler updates; 1 activation | 28,683.09–31,293.09 s |

The largest spread in final SMA increase is 0.23052 m, about 0.386% of the
mean increase. Thus these runs show closely similar SMA behavior, but not
numerically identical trajectories or link timing.

The v1 scheduler's activation counter reports two activations. Its schedule
log shows one extra initialization-time event: helper 2 is selected at t=0,
then deselected about 5.9 microseconds later. This is a near-zero SMA-score
roundoff effect for the helper that begins directly radial to the target; it
is not the physical link event listed in the table. After initialization, v1
also selects helper slot 20.

The prototype does not expose accepted-step link history. Its CSV therefore
re-evaluates the prototype's actual candidate selector at each saved solution
state; its interval and event counts are sampled diagnostics, not a count of
all RHS-level scheduler evaluations.

## Same and different

- **Same in this scenario:** all candidate links have identical force
  parameters, there is one target and one target terminal, and only one helper
  is in range at the useful encounter. Each implementation chooses the
  prograde, SMA-increasing helper (slot 20) and turns the link off once that
  positive opportunity ends.
- **Different scheduling cadence:** the prototype calls its SMA selection
  during RHS force evaluation. The v1 and 9 schedulers update at accepted
  solver steps. This moves their activation/deactivation times by several
  seconds and changes the integrated SMA by fractions of a metre.
- **Different optimizer scope:** the prototype and v1 select one helper for
  one target. The 9 implementation scores eligible graph edges and solves a
  maximum-total-score terminal matching. In this single-target star graph,
  the matching reduces to choosing one helper.
- **Different score scaling:** prototype and v1 rank helpers using the
  direction-dependent GVE score, assuming equal force magnitudes. The 9
  implementation scores actual force and target mass. If candidate links had
  different power or magnification, the rankings could differ.
- **Different availability checks:** the prototype checks line of sight and
  range; v1 checks range; 9 requires positive separation and separation no
  greater than the registered edge's own maximum range. Positive separation
  excludes coincident spacecraft (zero distance); it is not an additional
  minimum-range setting. Every link in this run used the same 200 km maximum.
  The model also checks its forbidden-pair set, but this run constructed the
  model with the default empty set, so that check excluded no links. The
  simulation passed spacecraft activity flags, but all spacecraft were
  active; the activity check likewise excluded no links. No candidate was
  occulted in this geometry, so the prototype's extra line-of-sight check
  excluded no links either.
- **Small source-constant difference:** the v1 GVE scoring routine uses its
  source constant 3.986004418e14 m^3/s^2, while the propagation in this
  comparison uses 3.98600436233e14 m^3/s^2. This does not change the selected
  helper here, but means the implementations are not perfectly parameter
  identical.

The 4/v1 project snapshot does not bind its exported laser scheduler into the
`SpaceAGORA` package namespace. Its runner therefore loads and executes the
actual `laser_link_effectors.jl` source in a small wrapper and supplies the
same central-gravity ODE. The prototype runner uses its existing
`run_open_cavity_multi`; 9 uses its native `run_simulation` path.

## Reproduce

Run these commands from the workspace root:

```sh
julia --project="9_SpaceAGORA.jl-main_full integration" SciTech_Paper/run_9_full.jl
julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/run_4_v1.jl
julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/run_1_prototype.jl
```

The per-run `summary.toml` and `schedule.csv` files are in `results/9_full`,
`results/4_v1`, and `results/prototype`. The prototype runner gives the legacy
runner its own `runs/prototype/plot_cache` directory because that runner
clears its supplied image directory at startup, even when plot generation is
disabled.

All three runners also accept two positive integer arguments: helper count
and initial target orbital-period count. For example, append `50 50` to run
50 helpers for 50 periods. Explicit arguments save results separately under
`results/50_helpers_50_orbits/` (and run scratch data under the corresponding
`runs/50_helpers_50_orbits/` directory), preserving this default comparison.

The completed larger run is documented in
[50-helper, 50-orbit comparison](SMA_SCHEDULER_COMPARISON_50_HELPERS_50_ORBITS.md).

All three runners additionally accept `schedule target_inclination_deg`
after those two counts. For example, `50 50 gve_inc 1` runs the inclination
scheduler with an initially one-degree target orbit. These cases use a
separate output directory including the policy and inclination. See the
[inclination experiment](INCLINATION_COMPARISON_50_HELPERS_50_ORBITS.md).
The 9 runner uses the newly added physical-rate `:gve_inc` policy.

## Code used

- [9 integration runner](./run_9_full.jl)
- [4 / v1 runner](./run_4_v1.jl)
- [Prototype runner](./run_1_prototype.jl)
- [9 SMA scoring and matching](../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/2_dynamic_programming_models.jl)
- [9 graph availability](../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/1_constellation_graph_models.jl)
- [v1 laser scheduler and forces](../4_SpaceAGORA.jl-main_v1/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl)
- [Prototype laser scoring and selection](../1_Kuang's%20Prototype%20Code/functions/2_Laser_Forces_ver2.jl)
