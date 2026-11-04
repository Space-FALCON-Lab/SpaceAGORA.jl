# Range-only versus positive along-track: run settings

## Requested comparison

One target at **1000 km**, **200 helpers at 1050 km**, **0.5-degree initial
relative inclination**, and **80 initial target orbital periods**, with **J2**.
Both policies start from identical states and use identical physical and
numerical parameters.

The unspecified initial conditions and laser parameters follow the existing
SciTech study defaults below. Confirmation was unavailable, so these choices
are recorded explicitly rather than treated as user-specified values.

## Initial conditions and physics

| Setting | Value |
|---|---|
| Target spacecraft index | 1 |
| Helper spacecraft indices | 2 through 201 |
| Target altitude above Earth's equatorial radius | 1000 km |
| Helper altitude above Earth's equatorial radius | 1050 km |
| Target initial inclination | 0.5 deg |
| Helper initial inclination | 0 deg |
| Initial eccentricity, all spacecraft | 0 |
| Initial RAAN, all spacecraft | 0 deg (conventional for equatorial helpers) |
| Initial argument of periapsis | 0 deg (conventional for circular initial orbits) |
| Target initial true-anomaly parameter | 1 deg |
| Helper phase rule | `360*(spacecraft_index-2)/200` deg |
| Mass of each spacecraft | 227 kg |
| Maximum laser range | 200,000 m |
| Laser input power | 10,000 W |
| Amplification factor | 100 |
| Transmission factor beta / efficiency eta | 1 / 1 |
| Speed of light | 299,792,458 m/s |
| Force on each endpoint of an active link | `eta*beta*amplification*power/c` N, approximately 0.003335640952 N |
| Earth gravitational parameter | 3.98600436233e14 m^3/s^2 |
| Earth equatorial radius | 6,378,136.6 m |
| J2 | 0.00108263 |
| Other forces | No drag, SRP, third-body gravity, or helper-helper lasers |
| Frame | Fixed Earth-equatorial inertial frame, J2 symmetry axis along +z |
| Epoch label | 2026-01-01T00:00:00; no epoch-dependent force models are used |

Every helper is dynamically propagated with central gravity and J2. The active
helper also receives the equal-and-opposite laser force; helpers are **not**
prescribed, unaffected trajectories. Only one target-helper link may be selected
at a time. Force direction on the target is helper-to-target.

The analytic J2 acceleration is a local, unchanged snapshot of the SciTech
gravity calculation, saved as [`Code/gravity_kernel.jl`](../Code/gravity_kernel.jl).
The simulation no longer includes source files from the shared library.
The sibling `9_SpaceAGORA.jl-main_full integration` project is used only as the
existing Julia dependency environment; its source and project files are not
edited. All comparison implementation and test-runner code is under `Code`.

## Policies

Implementation: [`scheduling_policies.jl`](../Code/scheduling_policies.jl).
Reference: `ORACLE_SpaceAGORA_Integration/2_SpaceAGORA.jl/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl`.

### Range-only: `naive_next_entering`

- At initialization, select the closest in-range helper.
- Keep that helper while it remains in range, even if another helper becomes closer.
- When the current helper leaves range, consider **only helpers newly entering
  range at that scheduler update**, choosing the closest.
- When idle, select the closest available helper only if the previous update had
  no in-range helpers; otherwise consider only newly entering helpers.
- Consequently, the policy can remain idle even while an already-in-range helper
  is available. This is intentional reference behavior, not a nearest-helper
  fallback.
- No along-track sign or orbital-element score is used.

### Positive along-track: `positive_along_track`

- Define the target's instantaneous RTN transverse direction as
  `T_hat = normalize(r cross v) cross normalize(r)`.
- A helper qualifies when it is in range and
  `dot(normalize(r_target-r_helper), T_hat) > 0`.
- Keep the current helper while these conditions remain true.
- When a replacement is needed, select the in-range helper with the largest
  **positive directional projection**. If none qualifies, remain idle.
- This is not global reoptimization at every step, velocity-direction scoring,
  or GVE/SMA-rate maximization.

Both policies use inclusive `distance <= range` eligibility, deterministic
helper-order ties, and scheduler updates at initialization and every accepted
solver step. There is no additional Earth-occultation gate, consistent with
the reference range flags. At these altitudes, a 200 km link cannot intersect
Earth.

## Numerical implementation and interpretation

| Setting | Value |
|---|---|
| Integrator | OrdinaryDiffEq Tsit5, adaptive |
| Relative / absolute tolerance | 1e-12 / 1e-12, SI state components |
| Maximum accepted time step | 1 s |
| Duration definition | `80*T0`, where `T0 = 2*pi*sqrt((R_earth+1,000,000)^3/mu)` |
| Initial target period T0 | 6,307.11893783881 s |
| Simulated duration | 504,569.51502710476 s = 5.83992494244 days |
| Saved trajectory grid | 200 samples per initial target period; 16,001 rows including both endpoints |
| Orbit-count axis | Elapsed time divided by the **initial** Keplerian target period |
| Orbital elements | Instantaneous osculating elements derived from Cartesian state and central mu |
| RAAN plotting | Unwrapped over every saved sample, in degrees |

**Force-integration distinction:** the two policy state machines reproduce
SpaceAGORA's reference behavior, but this comparison integrates laser forces
continuously in the ODE, rather than using the older runner's discrete velocity
kicks. Separation is checked at every RHS evaluation, so force is zero outside
range or at zero separation. The along-track sign is tested at scheduler updates,
not continuously between them. No promise of bit-for-bit trajectory equality
with the old runner is made.

With J2, the plotted elements include short-period variations. The orbit axis
is a fixed reference-period count, **not** an accumulated count of actual
revolutions or a time-varying osculating-period integral. Both plots display the
same samples with different horizontal units.

## Files and reproduction

The output directory hierarchy matches the `Simulation_absolute_RN` convention:

```text
ORACLE_SciTech_2027/
  Output_3_rangeOnly_vs_positiveAlongTrack/
    Code/
    Output/
      RUN_SETTINGS.md
      comparison_summary.csv
      images/orbital_elements/
      paper_plot_mode/
        h1050km_t1000km/
          ih0.0deg_it0.5deg/
            N200/
              T504570s/
                naive_next_entering_e0.0000_nu1.0000/
                positive_along_track_e0.0000_nu1.0000/
  Tests/
    Output_3_rangeOnly_vs_positiveAlongTrack/
      validation/
```

`T504570s` is the rounded duration used in the folder name; the actual
integration still ends at 504569.51502710476 s. Each policy leaf contains the
original comparison's Cartesian Feather, target-element CSV, schedule CSV,
and run metadata. This matches the **folder structure**, not the additional
absolute-RN campaign's encounter/RTN/work diagnostic schema: those diagnostics
were not recorded by this comparison and are not replaced with dummy values.
The target remains spacecraft 1, as recorded in the metadata.

- [`simulation.jl`](../Code/simulation.jl): settings, coupled propagation, orbital-element
  conversion, full-state output, and output verification.
- [`scheduling_policies.jl`](../Code/scheduling_policies.jl): both rebuilt policies and
  range-gated pair force.
- [`gravity_kernel.jl`](../Code/gravity_kernel.jl): local central/J2 acceleration.
- [`layout.jl`](../Code/layout.jl): shared code, output, and validation locations.
- [`plot_results.jl`](../Code/plot_results.jl): reloads verified CSV results and creates
  four-panel, standalone SMA, and SMA-change comparisons versus time and reference
  orbit count, all as PNG and PDF.
- [`test_policies.jl`](../Code/test_policies.jl): exact transition fixtures, 1,000-step
  comparisons against each original policy, geometry, force balance,
  conservation, J2 precession, step refinement, and output tests.
- [`verify_results.jl`](../Code/verify_results.jl): validates every saved production
  target sample, all spacecraft initial conditions, checksums and scheduler
  decisions, then compares both full-duration runs against 0.5 s maximum steps.
- `trajectory.feather` in each policy leaf: time, reference orbits, and all 201
  spacecraft's Cartesian positions and velocities, in SI units.
- `target_orbital_elements.csv` in each policy leaf: target elements and changes from
  the initial state, with time in seconds/days and initial-period orbit count.
- `schedule.csv` in each policy leaf: every scheduler decision, selected helper,
  distance, transverse projection, and number of in-range helpers. Blank distance
  and projection fields indicate no selected helper.
- `run_info.toml` in each policy leaf: exact settings, source/environment hashes,
  solver statistics, endpoint elements, output checksums, and completion status.
- `images/orbital_elements/orbital_elements_vs_time.{png,pdf}` and
  `images/orbital_elements/orbital_elements_vs_orbits.{png,pdf}`: target a, e, i,
  and RAAN, both policies overlaid.
- `images/orbital_elements/sma_vs_time.{png,pdf}` and
  `images/orbital_elements/sma_vs_orbits.{png,pdf}`: standalone semimajor-axis
  comparisons, both policies overlaid as solid lines in different colors, using
  the same saved samples.
- `images/orbital_elements/delta_sma_vs_time.{png,pdf}` and
  `images/orbital_elements/delta_sma_vs_orbits.{png,pdf}`: SMA-change histories
  `a(t) - a(0)` in km, subtracting each policy's first saved osculating SMA.
  Both curves start at zero and use solid lines. These are changes from the
  initial state, not differences against a laser-off control.
- `comparison_summary.csv`: final elements and changes for both policies.
- [`Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/`](../../Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/):
  test output, numerical step-refinement evidence, and source-isolation checks.

Run from `Output_3_rangeOnly_vs_positiveAlongTrack/` (the parent of this `Output`
directory) using Julia 1.12 and the existing sibling SpaceAGORA environment:

```bash
julia --startup-file=no --threads=1 Code/test_policies.jl
julia --startup-file=no --threads=1 Code/simulation.jl
julia --startup-file=no --threads=1 Code/verify_results.jl
GKSwstype=100 julia --startup-file=no --threads=1 Code/plot_results.jl
```

The scripts activate the sibling project automatically. No dependency upgrades
are requested. Simulation refuses to overwrite existing result directories.
To run another comparison, choose a new output root explicitly:

```julia
import Pkg
Pkg.activate("../9_SpaceAGORA.jl-main_full integration")
include("Code/simulation.jl")
using .ComparisonSimulation
run_comparison(settings=Settings(), root="/absolute/path/to/new_comparison/Output")
```

To regenerate figures from saved results, rerun only the plotting command.
The initial 10 s runs passed short two-orbit checks but the full-duration
10 s versus 5 s check found up to 16.66 m semimajor-axis difference and
6.49e-6 eccentricity difference for the range-only policy. Those coarse
results are retained under the dedicated Tests validation directory's
`initial_10s/` subfolder, not used in the final plots. Both production policies
were therefore rerun with 1 s maximum steps
and checked over the full 80 orbits against 0.5 s maximum steps. This is an
empirical refinement check, not a rigorous absolute error bound.

## Execution results

**Completed and verified on 2026-11-04 using Julia 1.12.1, one Julia thread.**
Both production integrations returned `Success`, reached exactly 80 reference
orbits, and saved 16,001 full-constellation samples. The machine-readable
`run_info.toml` files are the authoritative per-run records.

| Target quantity at 80 reference orbits | Range-only | Positive along-track |
|---|---:|---:|
| Semimajor axis, km | 7378.203082580 | 7386.880591111 |
| Change in semimajor axis, km | +0.066482580 | +8.743991111 |
| Eccentricity | 0.000366076567 | 0.000172253518 |
| Inclination, deg | 0.499686827757 | 0.498620082075 |
| Change in inclination, deg | -0.000313172243 | -0.001379917925 |
| Unwrapped RAAN, deg | -33.6170770848 | -33.4249237679 |

These are osculating endpoint values, including J2 variations, not
orbit-averaged values or laser-only differences against an unforced control.

Policy/physics/output checks and production/full-duration refinement checks
passed; current counts are recorded in the validation logs. The two four-panel
PNGs were visually inspected; PDF versions are also saved.

Maximum differences between the final **1 s** production runs and independent
**0.5 s** runs, evaluated at all 16,001 shared target sample times:

| Quantity | Range-only | Positive along-track |
|---|---:|---:|
| Semimajor axis, m | 3.825361 | 0.002192 |
| Eccentricity | 1.468158e-6 | 1.153423e-8 |
| Inclination, deg | 1.252738e-6 | 1.210993e-7 |
| Unwrapped RAAN, deg | 0.000525463 | 0.000294321 |

The acceptance thresholds were 10 m, 2e-6, 1e-4 deg, and 0.01 deg respectively.
The small range-only net semimajor-axis change is more sensitive to scheduling
time resolution than the much larger positive-along-track change; retain the
refinement figures above when interpreting that small endpoint difference.
Detailed evidence is in
[`full_duration_step_refinement.csv`](../../Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/full_duration_step_refinement.csv)
and [`production_verification.toml`](../../Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/production_verification.toml).

## Folder reorganization verification

The relocated runners passed **8,873 policy/physics/layout checks** and
**132,901 production/refinement checks**. Both 80-orbit cases were reproduced
under `Output/paper_plot_mode`. For each policy, `trajectory.feather`,
`target_orbital_elements.csv`, and `schedule.csv` are **byte-for-byte identical**
to the pre-reorganization results. Redundant old copies were removed only
after those comparisons passed; the original run metadata was retained under
the Tests validation directory's `prior_layout_metadata/` subfolder.

Checksums also confirmed that the shared library's source, extensions,
`Project.toml`, and `Manifest.toml` were unchanged during this reorganization.
Evidence:
[`relocation_data_check.log`](../../Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/relocation_data_check.log)
and [`shared_library_check.log`](../../Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/shared_library_check.log).
