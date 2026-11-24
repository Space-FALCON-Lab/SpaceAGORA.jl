# Range-based absolute R/N study (Output 8)

The default grid runs **only `:range_based`**. At initialization and after every
accepted solver step, the target (spacecraft `N+1`) selects the nearest eligible
helper within 200 km, without checking force direction. If none is in range,
no link is selected. Selection is held until the next accepted step.

The grid otherwise matches the [absolute R/N study](../../Output_2_Simulation_absolute_RN/Code/README.md):
five target altitudes (900, 950, 1000, 1050, 1100 km), three inclinations
(0, 0.5, 1 degree), and nine helper counts (1, then 100 through 800 in steps
of 100). This gives **135 cases**, each lasting **100 initial target periods**.
Helpers start at 1000 km in the equatorial plane. Laser power is 10 kW.

## Launch and monitor

From `ORACLE_SciTech_2027/Output_8_Simulation_absolute_RN_range_based`:

```bash
bash Code/launch_study.sh start           # detached tmux session
bash Code/launch_study.sh status
tail -f Output/simulation_progress_range_based.log
```

Use `run` instead of `start` for a foreground run. `stop-after-case` pauses
after the current case is verified; `resume` rechecks existing results and
continues. The batch also pauses if there is insufficient free space.

Results are saved under `Output/paper_plot_mode/` in case directories ending
in `range_based_e0.0000_nu1.0000`. Each contains the trajectory (1001 rows,
including signed and positive/negative target R/T/N delta-v), encounters, laser
firing intervals, range audit, longest encounter track, run metadata, and
verification report. Batch status and per-case logs are in `Output/_batch_range_based/`.
The copied SMA/INC comparison plots and paper-report scripts are not part of
this run. For encounter statistics from a completed grid, use
`julia Code/evaluate_encounters.jl`.

## Test without starting the grid

```bash
julia --startup-file=no --threads=1 ../Tests/absolute_RN/test_range_based_campaign.jl
```

This checks the default grid and runs a short test case through the production
worker and output verifier in a temporary directory under `Tests/`, not `Output/`.
