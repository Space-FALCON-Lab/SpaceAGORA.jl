# Last-two-orbit averages for Outputs 3, 4 and 5

These are additional versions of the existing orbital-element and RTN figures.
The original endpoint figures, source trajectories and endpoint summaries are
not replaced. No simulations are rerun and no force or scheduling settings change.

## Where the new figures are

Each folder below contains `80orbits/` and `100orbits/` subfolders:

- [Output 3: laser-on averages](Output_3_rangeOnly_vs_positiveAlongTrack/Output/images/last_two_orbits_average/)
- [Output 4: no-laser averages](Output_4_no_laser/Output/images/last_two_orbits_average/)
- [Output 5: laser-minus-no-laser averages](Output_5_Enable_Disable_Laser_compare/Output/images/last_two_orbits_average/)

Files start with `mean_last2orbits_`. Each plot has PNG and PDF versions,
a full-precision source CSV, and an explicit averaging-window title.
The existing colors, policy legends, helper-count axis, altitude columns,
inclinations, units and initial continuous-access thresholds are retained.
Axis limits are fitted to the averaged values, not copied from endpoint plots.
The averaged canvases are slightly wider to accommodate the no-laser panels'
longer tick labels without clipping; endpoint rendering dimensions are unchanged.

Outputs 3 and 4 each have orbital-element and laser-only RTN figures.
Output 5 also has the existing **actual velocity-state difference** figure,
which is distinct from the laser-only RTN integrals. There are 14 figures,
or **28 PNG/PDF files**, across the three outputs and two windows.

Each output's `Output/last_two_orbits_average_summary.csv` contains all its
averaged quantities. It is separate from `production_endpoint_summary.csv`.
Each new figure directory has a `figure_manifest.toml` recording the method,
window, case count, sample count and figure/source-CSV checksums.

## What "last two orbits" means

To keep these plots directly comparable with the existing 80/100-orbit figures,
an orbit is the same **initial Keplerian target period** used by the campaign:

\[
T_0 = 2\pi\sqrt{\frac{(R_E+h_{\rm target})^3}{\mu}}.
\]

| Figure window | Averaging interval | Inclusive source rows |
|---|---|---|
| 80 orbits | \(78T_0 \le t \le 80T_0\) | 781 through 801 |
| 100 orbits | \(98T_0 \le t \le 100T_0\) | 981 through 1001 |

The period, and hence the physical duration in seconds, depends on the initial
target altitude. This is **not** a measurement of the final two actual
revolutions, nor an integral of a changing osculating orbital period.
J2 and laser forces can change the phase relative to this reference clock.

The recorded runs contain 1001 uniform samples over 100 reference orbits:
10 intervals per orbit. Each averaging window therefore contains **21 samples
and 20 time intervals**, including both boundaries. Boundaries must match
saved times; the code does not choose a nearest endpoint or extrapolate.

## Averaging formula

For the already-computed plotted quantity \(q(t)\), the new point is

\[
\overline q =
\frac{1}{t_b-t_a}\int_{t_a}^{t_b}q(t)\,dt
\;\approx\;
\frac{\sum_{k=0}^{19}(t_{k+1}-t_k)
       [q(t_k)+q(t_{k+1})]/2}{t_{20}-t_0}.
\]

This is a **time-weighted trapezoidal mean**, not an equally weighted
arithmetic mean of all 21 samples. For the uniform saved grid it simplifies to

\[
\overline q =
\frac{\tfrac12q_0+q_1+\cdots+q_{19}+\tfrac12q_{20}}{20}.
\]

The two boundary samples receive half the weight of an interior sample.
Actual saved time differences are used in the implementation.

### Outputs 3 and 4: changes relative to initialization

- Orbital elements: average \(OE(t)-OE(0)\).
- Laser RTN: average each existing **cumulative** laser-only diagnostic
  \(\Delta v_R(t),\Delta v_T(t),\Delta v_{N+}(t),\Delta v_{N-}(t)\).
- Additional net, positive, negative and absolute R/T/N diagnostics are kept
  in the summary CSV even where the original figure shows only selected fields.
- Negative branches keep their signs.

In particular, the RTN average is **not** the additional impulse delivered only
within the last two orbits, \(\Delta v(t_b)-\Delta v(t_a)\).
It is also not an average acceleration or the spacecraft's full velocity.
The initial orbital element remains the reference; it is not replaced by the
orbital element at the start of the averaging window.

Output 4 has exactly zero laser RTN diagnostics. Its osculating orbital elements
can still differ from their initial values because J2 remains enabled.

### Output 5: differences between matched laser-on and no-laser runs

The instantaneous differences are recomputed from verified Output 3/4 source
trajectories, using the same definitions as the existing Output 5 products,
and then averaged:

- Orbital elements: \(\overline{OE_{\rm on}(t)-OE_{\rm off}(t)}\).
- Laser RTN integrals:
  \(\overline{\Delta v_{{\rm laser,on}}(t)-\Delta v_{{\rm laser,off}}(t)}\).
  The second term is exactly zero.
- Actual velocity difference: subtract the ECI velocity vectors at each
  matching time, project the result onto the **no-laser target's instantaneous
  RTN frame**, and then average each component.
- The norm panel is the **mean of the instantaneous difference magnitude**,
  \(\overline{\lVert\mathbf v_{\rm on}-\mathbf v_{\rm off}\rVert}\), not the
  magnitude of an averaged vector.

No pre-existing Output 5 difference trajectory is overwritten. The averaged
CSV records both source trajectories, their hashes and run IDs.

### Angles and nonlinear operations

Orbital elements are derived from each saved Cartesian state **before**
averaging. The code does not convert an averaged position/velocity into
orbital elements.

RAAN is unwrapped across the complete saved trajectory before subtracting
initial values or paired tracks and selecting the averaging interval.
Equatorial RAAN remains undefined (`NaN`) and is omitted from the RAAN panels;
it is never silently replaced by zero. Partially undefined or infinite values
are rejected.

## Metadata and units

The source CSV columns `orbits`, `time_s` and `row` identify the **end of the
averaging window**, for continuity with the endpoint schema. They do not imply
that the plotted value was taken at that single instant.

Additional columns explicitly record:

- `aggregation = time_weighted_mean` and `averaging_orbits = 2`;
- `window_start_orbits`, `window_end_orbits`;
- `window_start_row`, `window_end_row`;
- `window_start_time_s`, `window_end_time_s`;
- `sample_count = 21`.

CSV orbital-element units remain metres, dimensionless eccentricity, and
degrees. RTN/velocity units are m/s. Plot scaling is unchanged: semimajor axis
in km, eccentricity in \(10^{-4}\), inclination in \(10^{-3}\) degrees,
and RAAN in degrees.

## Reproduction and checks

From the `ORACLE_SciTech_2027` directory:

```bash
julia --startup-file=no --threads=1 \
  --project="9_SpaceAGORA.jl-main_full integration" \
  Output_3_rangeOnly_vs_positiveAlongTrack/Code/plot_last_two_orbits.jl
```

This one command generates averages for **all three outputs** and both windows.
It reads 162 laser-on runs and 81 no-laser runs. Existing verification,
trajectory checksums, exact time grids, initial states and paired settings
are checked before averaging. Missing or mismatched data cause an error.

Implementation:
[ComparisonResults.jl](Output_3_rangeOnly_vs_positiveAlongTrack/Code/ComparisonResults.jl)
and the shared
[paper renderer](Output_2_Simulation_absolute_RN/Output/images/for_paper/plot_paper_results.jl).

Tests:

```bash
julia --startup-file=no --threads=1 \
  --project="9_SpaceAGORA.jl-main_full integration" \
  Output_3_rangeOnly_vs_positiveAlongTrack/Code/test_last_two_orbits.jl unit

# After generating the figures:
julia --startup-file=no --threads=1 \
  --project="9_SpaceAGORA.jl-main_full integration" \
  Output_3_rangeOnly_vs_positiveAlongTrack/Code/test_last_two_orbits.jl saved
```

Tests cover nonuniform time weights, exact window boundaries, averages of
changes rather than endpoints, zero-laser diagnostics, paired differences,
undefined RAAN and invalid inputs. Saved-data checks independently form the
trapezoidal weights and compare every averaged quantity with the source
trajectories, including the existing Output 5 difference trajectories.
They also check source CSVs, checksums and every plotted data curve.

The generated figures passed **86 unit checks and 18,311 saved-data/plot checks**.
The original endpoint figures, figure-source CSVs, manifests and endpoint
summaries were also checksum-compared before and after generation and were
unchanged. The endpoint rendering path was separately regression-tested.

### Accuracy limitation

These averages use the existing 21 saved samples, not the solver's complete
accepted-step history or dense output. Trapezoidal integration is an
approximation to the continuous-time mean; short events or higher-frequency
variations between saved samples are not resolved. Two reference periods
reduce phase sensitivity but need not remove every J2 or laser-induced
oscillation. No continuous-time convergence claim is made.
