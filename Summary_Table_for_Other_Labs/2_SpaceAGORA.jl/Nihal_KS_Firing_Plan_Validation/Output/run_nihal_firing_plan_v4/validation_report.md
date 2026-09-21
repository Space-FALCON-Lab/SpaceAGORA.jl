# V4 Firing Plan Validation

Date: 2026-09-10

## Status

Schedule validation passed for all four supplied schedules. The simulation
pipeline loaded successfully in a fresh Julia process, and 25 focused tests
passed. All four controlled simulations and all four no-laser reference
simulations have now run using the uploaded initial conditions. Each trajectory
reached 31535.592545587482 seconds with 1001 saved samples. All 1356 post-run
output checks passed, including matching time grids, finite position/velocity
and orbital-element values, and required output-file presence.

## Final Target Changes

Values below are controlled minus no-laser reference at the final saved time,
not changes relative to the initial orbit. Angular differences are in degrees.

| Schedule | Delta a (km) | Delta e | Delta i (deg) | Delta RAAN (deg) | Delta omega (deg) |
| --- | ---: | ---: | ---: | ---: | ---: |
| max_a | +0.728489586 | +3.053805958e-6 | +0.000066623 | +0.000168300 | -0.006902901 |
| max_i | +0.090847911 | +2.861149541e-6 | +0.001738566 | +0.000073798 | +0.010319305 |
| max_joint | +0.723475617 | +2.797865085e-6 | +0.000065412 | +0.000163806 | +0.007182844 |
| max_omega | -0.004101592 | -2.348624469e-6 | +0.000001264 | +0.000014273 | +0.365044316 |

These checks establish run completion and output consistency, not agreement
with an independent Nihal/KS reference trajectory or optimizer predictions.

## Diagnostic Caveat

The saved `laser_active_helper_count` is zero only at the initial `t = 0`
sample for all four schedules and is one for all 1000 subsequent samples.
The existing engine registers the saving callback before the extra scheduler
callback, whose initializer activates the first scheduled helper. The initial
saved diagnostic therefore precedes scheduler initialization; it must not be
interpreted as an intentionally inactive first schedule interval. The post-run
checks explicitly verify this behavior. Shared engine code was not changed.

## Input Checks

| Schedule | Intervals | Spacecraft | Result |
| --- | ---: | ---: | --- |
| max_a | 1200 | 27 | Passed |
| max_i | 1467 | 27 | Passed |
| max_joint | 1209 | 27 | Passed |
| max_omega | 1211 | 27 | Passed |

All schedules begin at zero, have strictly increasing finite timestamps, and
contain binary helper flags with exactly one active helper per interval.
Target ID is 1; helper IDs are 2 through 27. Detailed interval measurements are
in [schedule_validation.csv](schedule_validation.csv).

The runner preserves every supplied switching time. The fixed mission end is
31535.592545587482 seconds, matching the previous v3 corrected/max_a manifest
and the last v4 start time plus the common 26.279660454656234-second base
interval. The extra switching times in three schedules must not shorten the
mission by using the average timestamp spacing, as the v2 loader would do.

The 25 tests covered exact timestamp and flag preservation, spacecraft mapping,
helper count and mission duration for every schedule, missing initial-condition
rejection, duplicate timestamp rejection, nonbinary flag rejection, and the
existence and row count of the generated validation CSV. Editor diagnostics
reported no errors in the changed Julia files.

## Initial Conditions Used

The run used the uploaded [initial_conditions.csv](../../Input/initial_conditions.csv),
passed explicitly to the runner. This resolved the missing-input blocker from
the first validation attempt. Initial conditions were not invented or inferred
from the schedules.

Required columns: `satellite`, `a_km`, `e`, `i_deg`, `raan_deg`, `omega_deg`,
`M_deg`. Each satellite ID from 1 through 27 must appear exactly once. Angles
are in degrees, semimajor axes are in kilometers, and `M_deg` is mean anomaly.

The runner's default remains `Input/run_nihal_firing_plan_v4/initial_conditions.csv`;
the actual shared `Input/initial_conditions.csv` location must be supplied as
the sole argument, as in the command below.

## Run Commands

From the workspace root, repeat schedule-only checks:

```bash
julia --startup-file=no --project=2_SpaceAGORA.jl \
  2_SpaceAGORA.jl/Nihal_KS_Firing_Plan_Validation/run_nihal_firing_plan_v4.jl \
  --check-inputs
```

Repeat the full run using the uploaded initial-condition file:

```bash
GKSwstype=100 julia --startup-file=no --project=2_SpaceAGORA.jl \
  2_SpaceAGORA.jl/Nihal_KS_Firing_Plan_Validation/run_nihal_firing_plan_v4.jl \
  2_SpaceAGORA.jl/Nihal_KS_Firing_Plan_Validation/Input/initial_conditions.csv
```

The runner reuses the existing v3 controlled/reference propagation and plotting
pipeline, including the 300 kg mass, J2 gravity, and external laser schedules.
Each schedule has written into its own `max_a`, `max_i`, `max_joint`, or
`max_omega` subfolder here, containing controlled CSV/Feather results, reference
CSV/Feather results, orbital-element time series, manifests, and 11 PNG plots
per schedule (44 total).