# Constellation terminal-gain scheduling

## Scope and provenance

This separate local repository integrates the fixed-reference terminal-gain
scheduler into SpaceAGORA. It is not the active-reference receding-horizon
scheduler, nor the battery, temperature, or data-constrained variant.

Base: `SpaceAGORA_MPC_Integration` at
`1af9cbab34452b277d9fd7efdd9038b36d3f4c9c`. On October 30, 2026, this exact
commit was confirmed directly against GitHub using
`git ls-remote upstream refs/heads/SpaceAGORA_MPC_Integration`.
The published integration branch is `Constellation-Controllability` in
`Space-FALCON-Lab/SpaceAGORA.jl`. The separate local integration checkout uses
the branch `constellation-controllability-interlinks`.
The older implementation is preserved on `constellation-controllability` at
`8ed59ffac`, based on `3ddcf3f`; it is not the current implementation.

The package name remains `SpaceAGORA`; no separate GitHub repository is needed.
In the local integration checkout, `upstream-local` points to the original SpaceAGORA checkout;
`upstream` points to `https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git`.

Scheduler source: the standalone fixed-reference scheduler archive, identified
by the following checksum for reproducibility.
SHA-256:
`8205f28b9d0cd604491785cd09fd473f55e51608ec3260ecdf3fa57296f953e3`.
The ZIP was extracted and run independently for regression comparison, not
included at runtime. This repository does not depend on that extraction or on
the original controllability repository.

## Run

From this repository, with Julia 1.12:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/unit/gnc/constellation_scheduling_tests.jl
julia --project=. test/unit/gnc/terminal_gain_interlink_engine_tests.jl
julia --project=. examples/Constellation_Terminal_Gain_Scheduling.jl --check
julia --project=. examples/Constellation_Terminal_Gain_Scheduling.jl --no-plots
# Plotting dependencies belong to the existing examples environment, not core.
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/Constellation_Terminal_Gain_Scheduling.jl
julia --project=. examples/Constellation_Interlink_Replay.jl
```

The full example runs four objectives, 27 satellites, five orbits and 1200
intervals. It writes `outputs/constellation_terminal_gain_scheduling/`, including initial
conditions, physical-time binary schedules, per-link terminal contributions,
linear and nonlinear element histories, terminal comparisons, and individual
schedule and element-comparison PDFs. The short check uses 12 intervals over
0.05 orbit and writes a separate `constellation_terminal_gain_scheduling_smoke/` folder.
Compilation is included in the first invocation, not a propagation benchmark.

## Code ownership and reuse

| File | Responsibility |
| --- | --- |
| Existing `src/models/dynamics/ks_dynamics/transform.jl` | Cartesian/KS conversion, kinematics and their differential maps |
| Existing `nonlinear_dynamics.jl` and `linearized_dynamics.jl` | KS point-mass/J2 dynamics and analytical KS state Jacobian |
| Existing `SimulationModel.FrameTransforms` | Orbital element conversion and RTN basis |
| Existing `src/models/vehicle/interlinks/` | Terminal graph, registration, range/occultation and laser forces |
| Existing `src/models/gnc/control/interlinks/interlink_scheduling.jl` | Candidate scoring and maximum-weight terminal matching |
| New `src/models/dynamics/ks_dynamics/physical_tangent.jl` | Acceleration input lift, physical-time tangent, constrained inverse differential and Cartesian push-forward |
| New `src/foundation/core/numerics/orbital_element_sensitivity.jl` | Analytical five-element output Jacobian |
| New `src/models/gnc/control/constellation_scheduling/constellation_scheduling.jl` | Reference/replay, ZOH gains and terminal-gain policy |
| New `src/orchestration/simulation/callbacks/terminal_gain_schedule_replay.jl` | Native engine configuration and exact interval-boundary stops |
| New `examples/Constellation_Terminal_Gain_Scheduling.jl` | Walker constellation configuration, CSV exports and plots |
| New `test/unit/gnc/constellation_scheduling_tests.jl` | Map, Jacobian, native geometry, ZOH and selection checks |
| New `test/unit/gnc/terminal_gain_interlink_engine_tests.jl` | Native engine emitter-switch/coast timing and KS trajectory agreement |

The existing drag-area input derivative is not a laser input matrix. The new
acceleration lift supplies the correct three inertial acceleration columns.
Likewise, the existing first-order Euler tangent is not substituted for the
scheduler's exact frozen-matrix ZOH. No second implementation of SpaceAGORA's KS
gravity/J2 dynamics or orbital/frame conversion is introduced.

The planner now constructs real `SpacecraftModel` objects, one virtual terminal
per spacecraft, an `InterlinkModel`, and native `LaserBeamLink` edges through
`register_link!`. All candidate edges share target terminal `(1,1)`, making
native terminal matching exactly the single-target at-most-one-link constraint.
It calls `update_availability!`, `interlink_wrench`, `score_candidates!` and
`select_interlinks!`; the nonlinear KS replay calls
`interlink_force_on_spacecraft`. No custom laser force law or duplicate
occultation implementation remains. `signed_link_schedule` is retained only
as the independent interval-argmax regression oracle, not the production solver.

`TerminalGainPolicy <: AbstractSchedulingPolicy` implements the native
`candidate_score` extension. Its scores are precomputed terminal gains, not
`SemimajorAxisRatePolicy` rates. The returned `plan.policy` masks out unselected
edges, so replay retains the planned activations rather than substituting a
different helper when a commanded link becomes unavailable. Native availability
and matching remain enforced. The geometry is evaluated continuously during
KS RHS calls; the native engine refreshes availability at accepted solver steps.

Native integration requires no changes to the existing engine or force effectors:

```julia
using SpaceAGORA
case = walker27_constellation()
times = collect(range(0.0, 5case.period_s; length=1201))
plan = plan_laser_schedule(case, times; objective=:a)
args = constellation_replay_configuration(case, plan.policy)
solution = run_simulation(args; return_solution=true,
    extra_callback_factories=(terminal_gain_schedule_callback,))
```

The callback factory uses `PresetTimeCallback` only to impose activation-boundary
stops. The existing interlink scheduler updates the selected graph at those times,
and the existing engine applies its laser effects. Before the policy starts and
at/after its final time, the policy returns zero. The native engine propagates
Cartesian states; the planner reference and its primary replay remain nonlinear KS.
The native replay example compares the two using identical schedules and each
route's own no-laser baseline. Its `--full` option runs five orbits; the default
is a short 0.05-orbit check.

Native undirected keys are canonical, for example `((1,1),(2,1))`, even though
the directed schedule labels that as emitter 2 illuminating target 1. The planner explicitly
selects the target endpoint's wrench; it does not assume the second key is the
receiver. CSV schedules retain the original emitter-to-target labels. All
spacecraft share the configured mass. If supplying custom initial states, keep
the native spacecraft initial conditions and the KS `initial_cartesian` columns
consistent. Classical derivatives require nonzero eccentricity and inclination.

## Units and force model

All internal lengths and velocities are metres and metres/second, times are
seconds and angles are radians. CSV display angles are explicitly marked in
degrees; the initial-condition CSV uses radians. Unlike the ZIP's kilometre
formulation with KS energy state equal to minus twice specific orbital energy,
SpaceAGORA uses

$$
y=[p^T,q^T,h,t]^T,\qquad \rho=p^Tp,\qquad
h=-\epsilon=\frac{\mu}{\|r\|}-\frac{v^Tv}{2},\qquad dt/ds=\rho.
$$

Therefore the acceleration input in fictitious time is

$$
B_{a,s}=\begin{bmatrix}0_{4\times3}\\
\tfrac12\rho\Lambda(p)^T\\-\rho v^T\\0_{1\times3}\end{bmatrix}.
$$

For edge $e=(s,1)$, define $d_e=r_1-r_s$ and $\hat d_e=d_e/\|d_e\|$.
The receiver acceleration is $u_e\zeta_e C_rP\hat d_e/(mc)$ and the
emitter acceleration is its negative. The nonlinear replay propagates all
27 spacecraft with these forces and J2. The same geometry gate is reevaluated
inside every RHS evaluation. It requires positive separation, range at most
200 km, and clearance of the entire finite line segment above Earth plus a
100 km opaque atmosphere. Planning uses an additional 1 km range margin and
checks 11 samples per interval. These samples are not a continuous-time
feasibility proof. Replay diagnostics use 20 interior/boundary samples per
interval and report observed violations rather than claiming guaranteed clearance.

The simulation defaults are $P=10$ kW, $C_r=100$, $m=300$ kg,
$\mu=3.986004418\cdot10^{14}$ SI, $R_E=6378136.3$ m and
$J_2=1.08262668\cdot10^{-3}$. The coefficient $C_r$ is retained from the
configured effective force law; no additional reflected thermal or battery model
is implied by it.

## Linearization and scheduling equations

The reference is one nonlinear, all-links-off KS constellation trajectory.
It is integrated with Tsit5, absolute/relative tolerances $2\cdot10^{-12}$,
maximum physical step one quarter scheduling interval, and exact stops at
schedule boundaries. The nonlinear controlled replay uses the same settings.
All clocks evolve with $dt/dt=1$ because this integration uses physical time,
not a separate fictitious-time endpoint solve.

At each interval midpoint, SpaceAGORA's analytical fictitious-time Jacobian
is converted to physical time:

$$
A_t=\frac{A_s}{\rho}-\frac{f_s}{\rho^2}\frac{\partial\rho}{\partial y},
\qquad B_t=B_{a,s}/\rho.
$$

With $H=\partial(r,v)/\partial y$ and the constrained lift $G$ satisfying
$HG=I_6$, the Cartesian matrices are

$$
A_C=(\dot H+HA_t)G,\qquad B_C=HB_t=\begin{bmatrix}0\\I_3\end{bmatrix}.
$$

The lift respects $h=-\epsilon$; its energy row is
$[-\mu r^T/\|r\|^3,-v^T]$. The RTN basis $M=[e_R,e_T,e_N]$ maps RTN
components to inertial components. With $R=\operatorname{blkdiag}(M^T,M^T)$,

$$
A_R=\dot R R^T+R A_C R^T,\qquad
B_{e,R}=R B_C\frac{C_rP}{mc}\hat d_e.
$$

The state components here are rotated inertial position and velocity
perturbations. They are not position and its derivative in a rotating frame.
J2 contributes to the frame derivative as well as the state Jacobian.

Freeze $A_R,B_R$ over each interval and evaluate

$$
\exp\left(\Delta t_k\begin{bmatrix}A_{R,k}&B_{R,k}\\0&0\end{bmatrix}\right)
=\begin{bmatrix}F_k&V_k\\0&I\end{bmatrix},\qquad
\delta x_{k+1}=F_k\delta x_k+V_ku_k.
$$

For $z=[a,e,i,\Omega,\omega]^T$, the output map is the analytical
$C_k=(\partial z/\partial(r,v))\operatorname{blkdiag}(M_k,M_k)$.
The code differentiates angular momentum, node vector and eccentricity vector
analytically, then uses the derivative of `atan2` for argument of periapsis.
It does not numerically difference propagated element errors.

The terminal contribution of activation $e$ in interval $k$ is

$$
g_{e,k}=C_NF_{N-1}\cdots F_{k+1}V_{k,:,e}.
$$

Backward multiplication computes these contributions once. The original MILP is

$$
\max\sum_{k,e}(w^Tg_{e,k})u_{e,k},\qquad
u_{e,k}\in\{0,1\},\quad u_{e,k}\leq\zeta_{e,k},\quad\sum_eu_{e,k}\leq1.
$$

It separates exactly by interval: choose the available edge with the largest
strictly positive score; otherwise choose none. Native terminal matching solves
this star graph without Python or a MILP dependency. Its deterministic tie rule
may choose a different tied optimum from the independent first-index oracle;
objective sums must still agree. Individual objectives use their physical row.
For the joint objective,

$$
w^Tg=\frac13\left(\frac{\Delta a}{1000\ \mathrm{m}}+
\frac{\Delta i}{\pi/180\ \mathrm{rad}}+
\frac{\Delta\omega}{\pi/180\ \mathrm{rad}}\right).
$$

The archived native-unit roundoff cutoffs are retained: $10^{-9}$ m for $a$,
$10^{-16}$ for $e$, and $10^{-14}$ rad for angles. They only suppress tiny
terminal contributions, not propagated dynamics. Selected contributions sum
to the terminal linear prediction; the example checks this identity.

## Interpretation and limits

Nonlinear changes mean controlled minus no-laser elements at the same physical
time, with angular differences wrapped to $[-\pi,\pi]$. They are not changes
relative to initial elements; J2 alone can change osculating elements.

At the zero-input reference, laser state derivatives are multiplied by zero.
Consequently the target-only linear model is the target block of the full
first-order constellation model. Emitter recoil is included in nonlinear
replay, but its effect on the subsequently applied beam direction is absent
from this first-order zero-input model. This integration preserves that original
limitation, including possible close-approach prediction discrepancies. It does
not claim the original method is accurate for arbitrary long horizons.

The unit tests check $HG$, Cartesian Jacobian equivalence, acceleration input
mapping, finite-difference element derivatives, occultation, ZOH and selection.
Full numerical ZIP regression results are recorded alongside the example output.

## Completed validation

Both the untouched archive and this integration were run for all four objectives
on the full 1200-interval, five-orbit grid. All 4800 interval decisions agree
exactly, including the emitter IDs. No sampled nonlinear geometry violations
were found. The table below gives the new results; it does not imply the
original first-order predictor is exact.

| Objective | Linear delta a (m) | Nonlinear delta a (m) | Linear objective-element change | Nonlinear objective-element change |
| --- | ---: | ---: | ---: | ---: |
| Maximize a | 728.701188 | 728.582296 | 728.701188 m | 728.582296 m |
| Maximize i | 150.840930 | 90.118411 | 0.001733471568 deg | 0.001738960051 deg |
| Maximize argument of periapsis | -4.697779 | -5.339423 | 0.364836121303 deg | 0.365100765417 deg |
| Joint | 724.230456 | 723.078030 | See all five rows in CSV | See all five rows in CSV |

Over all output times and objectives, maximum integrated-versus-archive
linear-history differences are $1.41\cdot10^{-8}$ m in a,
$2.28\cdot10^{-16}$ in e, and $1.68\cdot10^{-12}$ degrees in angles.
Maximum nonlinear-history differences are 0.001630 m in a,
$2.47\cdot10^{-10}$ in e, and $1.21\cdot10^{-7}$ degrees in angles.
The SI energy/state scaling differs from the archive's kilometre scaling, so
equal adaptive solver tolerances do not enforce identical integration steps
or dense interpolation. Bit-for-bit nonlinear histories are not expected.
The explicit nonlinear regression limits, retained from the previous integration,
are 1 cm, $10^{-9}$ in e, and
$10^{-6}$ degrees. Linear limits are independently much tighter: $10^{-6}$ m,
$10^{-12}$ in e, and $10^{-10}$ degrees. Exact schedule equality is required.

The detailed measurements are in
`outputs/constellation_terminal_gain_scheduling/archive_regression.csv`; terminal predictions and
nonlinear results are in `terminal_comparison.csv` in the same directory.
To repeat the independent comparison after running the untouched archive:

```bash
julia --project=. test/validation/compare_scheduler_baseline.jl /path/to/archive/output/exact
```

The archive comparison has 48 assertions. The new branch's existing reusable KS
dynamics and energy-convention regression files pass all 48 assertions. Existing
interlink scoring tests pass 24 assertions, and existing laser-physics tests pass
186 assertions. The focused integration tests pass 37 assertions, covering
physical clock agreement, native endpoint sign conventions, equal/opposite
acceleration, exact star matching, forbidden links and no replay substitution.
The native-engine switch/coast test passes 75 assertions, including exact stops
at 10, 20 and 30 seconds and all 27 spacecraft agreeing with KS to within 1 mm
position and 1 micrometre/second velocity at the end of that test.
The entire unrelated atmosphere, aerobraking and visualization suite was not run.

### Native engine five-orbit replay

`julia --project=. examples/Constellation_Interlink_Replay.jl --full` runs the
semi-major-axis plan through `run_simulation`, using the native interlink engine,
plus its own no-laser baseline. Every recorded active link was checked against
the same fixed plan. The comparison is in
`outputs/constellation_native_interlinks_5orbits/terminal_comparison.csv`, and the actual
engine selections are in `native_executed_schedule.csv` beside it.

Both nonlinear implementations give approximately 728.582296 m in delta a.
Native Cartesian minus nonlinear KS terminal differences are:

| Element | Native minus KS |
| --- | ---: |
| a | $4.38746\cdot10^{-6}$ m |
| e | $-1.23184\cdot10^{-12}$ |
| i | $6.63464\cdot10^{-12}$ deg |
| RAAN | $-6.24375\cdot10^{-10}$ deg |
| Argument of periapsis | $5.86282\cdot10^{-9}$ deg |

These are integration/coordinate-route agreement figures, not linear-prediction
errors. The native-engine full-horizon comparison covers the a objective; the
four-objective scheduler/archive regression uses nonlinear KS replay with the native
SpaceAGORA link geometry and force functions.
