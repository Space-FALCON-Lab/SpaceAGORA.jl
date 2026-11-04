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

The analytic J2 acceleration is reused directly from the sibling SciTech
library's
[`central_acceleration.jl`](../9_SpaceAGORA.jl-main_full%20integration/src/models/environment/gravity/central_acceleration.jl).
The compact six-state-per-spacecraft propagator avoids modifying the existing
SciTech campaign or its GVE matching scheduler.

## Policies

Implementation: [`scheduling_policies.jl`](./scheduling_policies.jl).
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

- [`simulation.jl`](./simulation.jl): settings, coupled propagation, orbital-element
  conversion, full-state output, and output verification.
- [`scheduling_policies.jl`](./scheduling_policies.jl): both rebuilt policies and
  range-gated pair force.
- [`plot_results.jl`](./plot_results.jl): reloads verified CSV results and creates
  both four-panel comparisons as PNG and PDF.
- [`test_policies.jl`](./test_policies.jl): exact transition fixtures, 1,000-step
  comparisons against each original policy, geometry, force balance,
  conservation, J2 precession, step refinement, and output tests.
- [`verify_results.jl`](./verify_results.jl): validates every saved production
  target sample, all spacecraft initial conditions, checksums and scheduler
  decisions, then compares both full-duration runs against 0.5 s maximum steps.
- `results/<policy>/trajectory.feather`: time, reference orbits, and all 201
  spacecraft's Cartesian positions and velocities, in SI units.
- `results/<policy>/target_orbital_elements.csv`: target elements and changes from
  the initial state, with time in seconds/days and initial-period orbit count.
- `results/<policy>/schedule.csv`: every scheduler decision, selected helper,
  distance, transverse projection, and number of in-range helpers. Blank distance
  and projection fields indicate no selected helper.
- `results/<policy>/run_info.toml`: exact settings, source/environment hashes,
  solver statistics, endpoint elements, output checksums, and completion status.
- `plots/orbital_elements_vs_time.{png,pdf}` and
  `plots/orbital_elements_vs_orbits.{png,pdf}`: target a, e, i, and RAAN, both policies overlaid.
- `comparison_summary.csv`: final elements and changes for both policies.
- `validation/`: test output and numerical step-refinement evidence.

Run from this directory using Julia 1.12 and the existing sibling SpaceAGORA
project/environment:

```bash
julia --startup-file=no --threads=1 test_policies.jl
julia --startup-file=no --threads=1 simulation.jl
julia --startup-file=no --threads=1 verify_results.jl
GKSwstype=100 julia --startup-file=no --threads=1 plot_results.jl
```

The scripts activate the sibling project automatically. No dependency upgrades
are requested. Simulation refuses to overwrite existing result directories.
To run another comparison, choose a new output root explicitly:

```julia
include("simulation.jl")
using .ComparisonSimulation
run_comparison(settings=Settings(), root="/absolute/path/to/new_comparison")
```

To regenerate figures from saved results, rerun only the plotting command.
The initial 10 s runs passed short two-orbit checks but the full-duration
10 s versus 5 s check found up to 16.66 m semimajor-axis difference and
6.49e-6 eccentricity difference for the range-only policy. Those coarse
results are retained under `validation/initial_10s/`, not used in the final
plots. Both production policies were therefore rerun with 1 s maximum steps
and checked over the full 80 orbits against 0.5 s maximum steps. This is an
empirical refinement check, not a rigorous absolute error bound.

## Execution results

Execution status and measured endpoints will be recorded here after both
production runs and plot verification complete. The machine-readable
`run_info.toml` files are the authoritative per-run records.
