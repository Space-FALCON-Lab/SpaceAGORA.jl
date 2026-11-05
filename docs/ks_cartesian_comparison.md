# Native KS versus Cartesian comparison

This is the local SpaceAGORA port of the nonlinear comparison and original
fixed-reference scheduling comparison. It deliberately uses this repository's
KS implementation, not copied dynamics from the standalone comparison ZIP.
The historical affine/GVE/PPTX studies are not silently substituted for the
paper scheduler, and have not been ported by this driver.

## Run and review

The integration branch is `Constellation-Controllability` in `Space-FALCON-Lab/SpaceAGORA.jl`.
The benchmark has its own environment because DP5 is provided by
`OrdinaryDiffEqLowOrderRK` rather than this version of `OrdinaryDiffEq`.
The source dependency points at this checkout, not a registered SpaceAGORA release.

```bash
julia --project=examples/ks_cartesian_comparison -e 'using Pkg; Pkg.instantiate()'
julia --project=examples/ks_cartesian_comparison test/validation/ks_cartesian_comparison_tests.jl
julia --project=examples/ks_cartesian_comparison examples/ks_cartesian_comparison/comparison.jl short
python3 -m venv /tmp/agora_comparison_plot_env
/tmp/agora_comparison_plot_env/bin/pip install -r examples/ks_cartesian_comparison/requirements.txt
/tmp/agora_comparison_plot_env/bin/python examples/ks_cartesian_comparison/plot_results.py
```

Python plotting needs matplotlib and pandas. `short` runs 1 and 5 orbits for
the nonlinear cases and all four five-orbit original scheduling objectives.
`full` extends the nonlinear horizons to 1, 5, 50, 100 and 200 orbits; it does
not extend the scheduling horizon. `nonlinear` and `schedule` run either half.
Full runs overwrite the CSV files in the same output directory, so retain a
copy of short results before running them. There are no battery, temperature,
data constraints or receding-horizon acceptance tests in this fixed-reference comparison.

To rerun **only** the nonlinear comparison at all five horizons while preserving
the short-run data and five-orbit schedules:

```bash
julia --project=examples/ks_cartesian_comparison examples/ks_cartesian_comparison/comparison.jl extended
/tmp/agora_comparison_plot_env/bin/python examples/ks_cartesian_comparison/plot_results.py
```

This writes `outputs/ks_cartesian_comparison/extended_horizons/`. When present,
the plotting/reporting script requires its complete 240-configuration sweep
and uses it instead of the older short-run nonlinear timings. Combined CSVs
are exported as `nonlinear_all_horizons.csv` and
`truth_refinement_all_horizons.csv` in the main results directory.

Trajectory sampling now fills a preallocated 6-by-N-by-time array rather than
concatenating thousands of small arrays with variadic `cat`. The sampled
values are unchanged (checked for exact equality), but the previous assembly
had substantial allocation/GC overhead at long horizons. All five horizons
are remeasured with the same new implementation; old total timings must not
be spliced into the new sweep. The interrupted pre-change run is preserved in
`interrupted_preallocation_baseline/` locally and excluded from the reports
and published results.

Outputs are under `outputs/ks_cartesian_comparison/`: initial conditions,
nonlinear errors/timings, truth convergence, terminal element changes,
schedule CSVs, histories, timing breakdowns, and PDF/PNG comparisons.

## Nonlinear inertial propagation

The ten spacecraft are original Walker IDs 1,4,7,10,13,15,18,21,24,27.
Their original shell/plane/slot geometry is retained; eccentricity is set to
0.01 or 0.3, with $a_i(e)=a_i(0.01)\,0.99/(1-e)$ preserving perigee
distance, exactly as in the earlier comparison. Each case uses its own target
Kepler period. Only point-mass gravity and J2 are active in this part.
All lengths are SI. Cartesian dynamics are

$$\dot r=v,\qquad\dot v=-\mu r/\|r\|^3+a_{J2}(r).$$

`ks_j2_acceleration_si` supplies the same perturbation in both routes.
KS uses SpaceAGORA's state $y_i=[p_i,q_i,h_i,t_i]$, with $h_i=-\epsilon_i$,
and calls `nonlinear_ks_rhs` directly. It does not call any Jacobian routine.
Writing $\rho_i=p_i^Tp_i$, the shared target independent variable is

$$\frac{dy_i}{ds_1}=\frac{\rho_1}{\rho_i}f_{s_i}(y_i),
\qquad\frac{dt_i}{ds_1}=\rho_1.$$

The factor is required even without interlinks. It makes all spacecraft clocks
advance together despite their different instantaneous orbital radii.
A continuous terminal event stops the solve at $t_1=t_f$. Dense interpolation
and scalar bisection recover the common physical output epochs. The benchmark
includes this endpoint event in solve timing and dense sampling in total timing;
it does not hide clock costs or claim to eliminate endpoint root finding.

Solvers: Tsit5 and DP5. Adaptive absolute and relative tolerances are both
$10^{-8},10^{-10},10^{-12}$. Fixed Cartesian steps are 5,20,100 seconds;
native KS takes $\Delta s=\Delta t_{\rm nominal}/a_1$, which is **not** a fixed
physical-time step. Results record the setting and coordinate route explicitly.

Truth is Cartesian Vern9 at absolute/relative tolerance $10^{-14}$, compared
against a separate $2\cdot10^{-14}$ solve. These are finite-precision reference
trajectories, not exact truth; `truth_refinement.csv` records the measured
uncertainty. Output sampling is 240 points per orbit plus the initial point.
Maximum and RMS errors use the Euclidean position/velocity differences over
all spacecraft and epochs. Velocity errors are in m/s, not mm/s.

The driver warms each solver/route before measurement and sets BLAS to one
thread. Solve timings include problem construction and integrator execution.
Total timing additionally includes coordinate conversion, physical-clock
sampling and error calculation, but excludes truth generation, compilation,
CSV writing, plotting and environment startup. Allocated bytes and GC times
are recorded, not subtracted. One warmed measurement per setting is a local
validation, not a statistically robust speed ranking.

## Original scheduling linear model

`plan_laser_schedule` now accepts `linearization=:ks` (unchanged default) or
`:cartesian`. No scheduling objective, availability constraint, cut-off,
selection rule or default propagation setting is changed.

The direct Cartesian matrices are

$$A_C=\begin{bmatrix}0&I\\
\mu(3rr^T/\rho^5-I/\rho^3)+\partial a_{J2}/\partial r&0
\end{bmatrix},\qquad B_a=\begin{bmatrix}0\\I\end{bmatrix}.$$

They use `ks_j2_acceleration_jacobian_si`, SpaceAGORA's existing analytical J2
derivative, but do not obtain $A_C$ by a KS coordinate push-forward.
For the KS route, the existing `ks_physical_tangent` gives

$$A_t=A_s/\rho-f_s(\partial\rho/\partial y)/\rho^2,
\quad B_t=B_{a,s}/\rho,
\quad A_C=(\dot H+HA_t)G,\quad B_a=HB_t.$$

Both routes then use the existing RTN transform, the same frozen-matrix
augmented exponential for ZOH, and the same analytical orbital-element map:

$$\delta x_{k+1}=F_k\delta x_k+V_ku_k,
\qquad g_{e,k}=C_NF_{N-1}\cdots F_{k+1}V_{k,:,e}.$$

All four plans are computed independently using the native terminal matcher,
on a shared nonlinear no-laser KS reference. This isolates Jacobian construction
from reference-integration errors. This linear comparison is the original
physical-time RTN paper workflow, **not** native-s KS ZOH and not an inertial
affine approximation to the full orbit. It has 1200 intervals over five orbits,
27 satellites, 10 kW and reflection/magnification parameter 100.

For each objective the KS-selected binary command is replayed using both
nonlinear KS physical-time integration (existing `replay_laser_schedule`) and
direct Cartesian integration. Both call the existing native interlink force and
availability APIs, including endpoint recoil. They impose the same activation
boundaries. Element changes subtract each route's own no-laser trajectory;
angular differences are wrapped. If the independently optimized schedules
differ, that count is reported and their CSVs retained, while the nonlinear
cross-coordinate comparison still uses one identical commanded schedule.

The all-off reference has no laser-dependent state derivatives. Therefore the
original target-only first-order prediction still omits geometry feedback from
laser-induced emitter/target displacement. Agreement between coordinate routes
does not repair this modeling limitation or guarantee linear/nonlinear agreement
for the inclination-maximizing close approach.

## Native inertial linear comparison

Each selected schedule is also propagated with two additional target-only
inertial ZOH models, saved in `maximize_*/inertial_zoh.csv`. These are a
reference-plus-perturbation model, not an affine approximation of the full orbit.
No RTN conversion is used in these predictions. The native KS matrices are
`evaluate_ks_dynamics_and_jacobians(...).state` and
`ks_acceleration_input_matrix` multiplied by the reference laser acceleration
columns. Their zero-input reference is the same nonlinear trajectory used by
the optimizer. The Cartesian model uses the direct inertial matrices above.

KS freezes the matrices at the physical interval midpoint and obtains the
nominal fictitious-time increment by Simpson quadrature:

$$\Delta s_k\simeq\frac{\Delta t_k}{6}
\left(\frac{1}{\rho(t_k)}+\frac{4}{\rho(t_{k+1/2})}
+\frac{1}{\rho(t_{k+1})}\right).$$

The exponential is evaluated after a similarity scaling with diagonal scales
$[\sqrt{a_1}I_4,\sqrt{\mu}I_4,\mu/a_1,\sqrt{a_1^3/\mu}]$ to avoid SI-state
conditioning problems. This scaling does not change the continuous equations.
After advancing the perturbation in $s$, the clock perturbation is removed:

$$\delta y_{k+1,t}=\delta y_{k+1,s}
-\frac{f_s(\bar y_{k+1})}{\rho_{k+1}}\delta t_{k+1,s},
\qquad\delta x_{k+1}=H_{k+1}\delta y_{k+1,t}.$$

This imposes the physical switching boundary and resets only the clock
perturbation, not the physical position or velocity perturbations. Equivalently,
the output correction is $H-f_C e_t^T$. Native-s and Cartesian frozen-matrix
discretizations need not agree exactly at finite steps. The native prediction
does not re-optimize the schedule and retains the original zero-input model's
geometry-feedback limitation. Here only the target tangent is propagated;
helpers enter the reference input geometry. It must not be described as a
coupled constellation active-reference STM.

`inertial_refinement.csv` repeats the prediction with 1, 4 and 16 linear
substeps per scheduling interval. The binary command and nonlinear reference
are unchanged. This separates frozen-matrix discretization error from the
original model's laser-geometry feedback error. The known selected input is
combined before the augmented exponential, which is algebraically identical
to forming all $V_e$ columns and then multiplying by the binary command.

## Verification boundaries

Tests check coordinate round trips, shared clocks, direct-versus-mapped
Jacobians, input lifts, Tsit5/DP5 agreement at common physical times, all four
short schedules, and nonlinear laser/recoil replay agreement. The existing
scheduling regression suite is also retained. No new comparison is presented
as a rerun of the historical 200-orbit experiments or matched-accuracy search.

## Initial short validation

The short run contains 96 nonlinear configurations, four truth-refinement
checks, four independently computed five-orbit schedule pairs, and inertial
ZOH refinements at 1, 4 and 16 substeps. The new suite passed 104 checks;
the pre-existing scheduling suite passed 37 checks. All 4800 decisions match
between the two scheduler routes and match the previously validated outputs.
The maximum common-clock endpoint mismatch was $1.46\cdot10^{-11}$ seconds.

For the semi-major-axis objective, the physical-time paper models both predict
728.701188299 m. Nonlinear KS gives 728.582295792 m; Cartesian gives
728.582300561 m. Across all four objectives the maximum difference in the
two linear delta-a results is $4.37\cdot10^{-11}$ m. The maximum difference
between their nonlinear delta-a results is 0.008192 m (inclination objective).
These are validation differences at the configured tolerances, not exact equality.

The inclination objective retains the original model discrepancy: linear
delta a is 150.840930 m, versus 90.118411 m for nonlinear KS and 90.126603 m
for nonlinear Cartesian. The port does not repair the absent laser-geometry
feedback in the zero-input first-order model.

The independent inertial comparison illustrates a different issue: frozen
coordinate matrices are not equivalent discrete models at finite steps.
For the semi-major-axis schedule, native-s KS / Cartesian inertial predictions
are 728.670125 / 669.576619 m at one substep, 728.692776 / 724.912554 m
at four, and 728.694191 / 728.457584 m at sixteen. The cross-coordinate gap
decreases by approximately sixteen each time the step is divided by four,
consistent with second-order midpoint freezing. Sixteen substeps is not a
claim of full convergence. The original RTN optimizer remains unchanged.

See `outputs/ks_cartesian_comparison/README.md` for tables of the fastest
passing tested nonlinear settings and every terminal element. Timings are
single warmed measurements: near-ties are not evidence of a reliable winner.

## Completed 50-200 orbit extension

The extended sweep is complete: 240 nonlinear configurations covering both
eccentricities, all five horizons, both solvers and both coordinate routes.
There are ten separate Vern9 reference-refinement checks. The updated native
comparison suite passes 112 checks, including exact agreement of old and new
sample-array assembly. All four recorded error metrics also agree exactly for
the 24 completed pre-change 50-orbit runs. Their median total allocation fell
from 2.503 GB to 0.193 GB without changing integration or sampled values.

Fastest **tested** settings satisfying maximum position error plus reference
refinement difference at most 0.01 m, and the analogous velocity criterion
at most 0.00001 m/s, gave these total runtimes (seconds):

| Eccentricity | Orbits | Cartesian | KS target s |
| --- | ---: | ---: | ---: |
| 0.01 | 50 | 0.3494 | 0.08981 |
| 0.01 | 100 | 0.7575 | 0.1873 |
| 0.01 | 200 | No tested setting passes | 1.0657 |
| 0.3 | 50 | 1.0268 | 0.09553 |
| 0.3 | 100 | 1.7428 | 0.2398 |
| 0.3 | 200 | No tested setting passes | 0.5441 |

These totals include sampling and metrics, not compilation, truth generation
or plotting. Solve-only times, GC times, allocated bytes, settings and all
errors are in `extended_horizons/nonlinear.csv`. This is a finite settings
grid with one warmed timing per configuration, not a universal speed ranking.
At 200 orbits, the smallest tested Cartesian position errors were 0.01564 m
(e=0.01) and 0.01871 m (e=0.3), both using 5-second steps. Smaller steps or
tighter tolerances have not been searched. No matched-accuracy speedup ratio
is asserted when a route has no passing tested setting.

The selected KS 200-orbit position errors were 0.0006492 m and 0.004042 m,
respectively. The largest reference-refinement difference over the entire
sweep was 0.0008515 m and 0.0000007563 m/s; the maximum common-clock endpoint
mismatch was 0.0000000009314 seconds. The measured refinement difference is
a safety margin, not a rigorous bound on the reference's true error.

Only the nonlinear horizons were extended. The original scheduling and
native inertial linear experiments remain the previously validated five-orbit
cases. The implementation, completed results and figures are included with
the comparison on the integration branch.
