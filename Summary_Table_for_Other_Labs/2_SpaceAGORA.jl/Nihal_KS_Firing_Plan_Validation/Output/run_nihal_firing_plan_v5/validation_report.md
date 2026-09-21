# V5 Firing Plan Validation

Date: 2026-09-11

## Results

The supplied v5 `max_joint` schedule was run with the shared
[initial_conditions.csv](../../Input/initial_conditions.csv). Both the controlled
trajectory and the no-laser reference completed with solver status `Success`.
Controlled runtime was 24.40 seconds; reference runtime was 18.54 seconds.
The full console output is in [run.log](run.log).

The schedule contains 1223 intervals for target 1 and helpers 2 through 27.
All flags are binary with exactly one helper active per interval. All switching
times were preserved, including the irregular intervals. The mission end is
31535.592545587482 seconds, using the same v4 duration: the final schedule start
plus the common 26.279660454656234-second base interval.

Both trajectories contain 1001 samples on identical time grids, spanning zero
through the mission end. All 340 post-run checks passed: sample counts, time
grids, finite position/velocity values for all spacecraft, finite orbital
elements, helper-count diagnostics, and expected nonempty output files.

Final target changes below are controlled minus no-laser reference, not changes
relative to the initial orbit:

| Quantity | Final Difference |
| --- | ---: |
| Semimajor axis | +0.397786004 km |
| Eccentricity | +2.158026757e-6 |
| Inclination | +0.000949776 deg |
| RAAN | +0.000070552 deg |
| Argument of periapsis | +0.186761327 deg |

These checks establish completion and output consistency, not independent
agreement with Nihal/KS trajectories or optimizer predictions.

## Output Files

- [Schedule validation](schedule_validation.csv)
- [Controlled CSV](max_joint/simulation_results.csv)
- [Controlled Feather](max_joint/simulation_results.feather)
- [Reference CSV](max_joint/reference/simulation_results.csv)
- [Reference Feather](max_joint/reference/simulation_results.feather)
- [Orbital-element time series](max_joint/oe_timeseries.csv)
- Eleven PNG plots in `max_joint/images`, plus controlled/reference manifests.

As in v4, the saved `laser_active_helper_count` is zero only at `t = 0` because
the saving callback precedes scheduler initialization; all 1000 later samples
report one active helper. This is a diagnostic initialization caveat, not an
intentionally inactive first schedule interval. The output checks explicitly
account for it. Shared simulation code was not changed.

## Reproduce

From the workspace root:

```bash
GKSwstype=100 julia --startup-file=no --project=2_SpaceAGORA.jl \
  2_SpaceAGORA.jl/Nihal_KS_Firing_Plan_Validation/run_nihal_firing_plan_v5.jl
```

Append `--check-inputs` for validation without propagation, or pass a different
initial-condition CSV path as the sole argument. The runner reuses the existing
v4 loader and v3 simulation pipeline, including 300 kg spacecraft mass, J2
gravity, and externally scheduled laser links. V4 outputs were not modified.