# Laser-enabled minus no-laser comparisons

Additional [last-two-orbit averaged figures](../Output/images/last_two_orbits_average/)
are available for the 80- and 100-orbit windows, including both laser RTN
integrals and actual velocity-state differences. See the shared
[averaging method and reproduction guide](../../LAST_TWO_ORBITS_AVERAGES.md).
The original endpoint figures and source simulations are unchanged.

[plot_results.jl](plot_results.jl) reads verified, checksum-matched trajectories
from Output 3 and Output 4. The combined Output 3 background job runs this step
automatically after both complete campaigns. It does not rerun propagation.

There are **162 pairs**, one per Output 3 case. Pairing requires identical:

- target/helper initial geometry, count, IDs and every initial Cartesian state;
- J2/physical parameters, numerical tolerances, cadence, duration and source hashes;
- all **1001 exact saved times**, with no interpolation or nearest-time matching.

Only scheduler, laser power and its derived force may differ. The baseline must
be the matching `no_laser`, zero-power case. Missing, incomplete, mismatched or
unverified cases stop processing rather than silently producing partial plots.

## Three distinct quantities

1. **Orbital elements:** osculating `OE_laser(t) - OE_no_laser(t)`. Because initial
   states are identical, this also equals the difference of their changes from
   initialization. RAAN is unwrapped separately over every saved row before
   subtraction. Equatorial RAAN is undefined, not assigned a false zero.
2. **Laser-only RTN integrals:** Output 2's signed R, signed T, positive N and
   negative N diagnostics, laser-on minus no-laser. The no-laser integrals are
   exactly zero, so these differences equal the laser-on integrals. Negative
   branches keep their signs. This is **not** a velocity-state difference.
3. **Actual target velocity difference:** subtract `v_laser(t)-v_no_laser(t)` in
   the shared inertial frame, then project that vector onto the **no-laser target's**
   instantaneous RTN basis. Do not subtract velocity components expressed in
   two different local frames. The fourth row is the ECI difference magnitude.
   This includes gravitational trajectory changes induced by the laser, not
   just its integrated acceleration.

## Results

`Output/paper_plot_mode/` mirrors each laser-on case's hierarchy. Every pair has
`laser_minus_no_laser.feather` with all 1001 target-difference samples and
`comparison_info.toml` with both source paths/run IDs/checksums, definitions,
matching checks and the resulting data checksum.

`Output/production_endpoint_summary.csv` preserves every quantity at 80 and 100
initial target periods with source provenance.

`Output/images/for_paper/` contains the 100-orbit figures, each as PNG and PDF,
with full-precision CSVs and a figure manifest:

- `laser_minus_no_laser_OE_from_feather_combined_e0.0000_nu1.0000_100orbits_h_cols`
- `laser_minus_no_laser_dv_RTN_split_N_from_feather_combined_e0.0000_nu1.0000_100orbits_h_cols`
- `velocity_difference_RTN_laser_minus_no_laser_100orbits_h_cols`

The same figures at 80 orbits are in `Output/images/80orbits/`. All retain three
altitude columns, four quantity rows and helper-count curves colored by initial
inclination and scheduler, like the original paper figures. Axis limits are
fitted to the new data; no old hard-coded limits hide negative or larger effects.
OE's dashed access thresholds remain an initial fixed-circular range-only estimate,
not a guarantee about a scheduler's actual firing coverage.

Regenerate from the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Output_5_Enable_Disable_Laser_compare/Code/plot_results.jl
```

Optional positional arguments are `LASER_OUTPUT_ROOT NO_LASER_OUTPUT_ROOT
DIFFERENCE_OUTPUT_ROOT`. Derived products can be regenerated; input simulation
files are never changed. Data handling/rendering is shared in
[ComparisonResults.jl](../../Output_3_rangeOnly_vs_positiveAlongTrack/Code/ComparisonResults.jl).
