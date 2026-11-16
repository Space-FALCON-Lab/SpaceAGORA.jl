# Matched no-laser production campaign

Additional [last-two-orbit averaged figures](../Output/images/last_two_orbits_average/)
are available for the 80- and 100-orbit windows. See the shared
[averaging method and reproduction guide](../../LAST_TWO_ORBITS_AVERAGES.md).
The original endpoint figures and source simulations are unchanged.

The [combined SSH-independent launcher](../../Output_3_rangeOnly_vs_positiveAlongTrack/Code/README.md)
runs this campaign automatically after Output 3. There are **81 cases**:
helpers at 1000 km/0 deg; targets at 1000/1050/1100 km and 0/0.5/1 deg;
helper counts 1/100/200/300/400/500/600/700/800; circular initial orbits,
target anomaly 1 deg; 100 initial target periods.

The shared Output 2 production engine, initial spacecraft ordering (`target=N+1`),
J2 gravity, tolerances, 10 s step cap, diagnostics and 1001-row output grid are
identical to the paired laser-on cases. The only changed physical setting is
laser power **0 W**; the `no_laser` scheduler selects no links. All helpers are
still propagated and all-pair encounters are still recorded.

Each no-laser case supplies both Output 3 schedulers with the same matched
baseline. Baselines are independently propagated for every helper count, not
copied from a different state-vector size.

Saved hierarchy:

```text
Output/paper_plot_mode/h1000km_t1050km/ih0.0deg_it0.5deg/N100/T...s/
  no_laser_e0.0000_nu1.0000/
    trajectory.feather
    encounters.feather
    laser_firing.feather
    firing_range_audit.feather
    longest_encounter.feather
    run_info.toml
    verification.toml
    progress.toml
```

Firing/audit tables are genuinely empty, with their normal schemas; integrated
laser impulse, work and RTN diagnostics are exactly zero. They are not fabricated
stand-ins for missing measurements. Cartesian trajectories and orbital elements
still evolve under J2. Verification checks these properties and conservation.

Standalone foreground reproduction, from the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Output_4_no_laser/Code/run_grid.jl
julia --startup-file=no --threads=1 --gcthreads=1 Output_4_no_laser/Code/plot_results.jl
```

For SSH-independent execution use the combined launcher, rather than these
foreground commands. Status is in `Output/_batch_no_laser/status.toml`; per-case
logs and a readable progress log use the same layout as Output 3.

Figures, full-precision endpoint/source CSVs and manifests go to
`Output/images/for_paper/` (100 orbits) and `Output/images/80orbits/`.
The RTN figure is correctly all zero; the OE figure shows total J2-only changes.
The same encounter evaluation and production case/endpoint summaries as
Output 3 are generated without overwriting any laser-on data.
