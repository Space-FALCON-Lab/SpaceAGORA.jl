# GRAM Trajectory Lookahead Fidelity Study

This study compares the GRAM trajectory-generation method against direct GRAM
point evaluation at every whole second over a fixed trajectory horizon while
sweeping the number of trajectory interpolation points and interpolation
methods. It writes one row per planet, method, and point count, plus averaged
error metrics and plots with separate planet/method traces for:

- mean wind-vector absolute error magnitude (m/s)
- mean density absolute error (kg/m^3)
- mean temperature absolute error (K)

By default, the study samples the generated path once per second across a
900-second trajectory horizon, compares each sample against direct GRAM point
evaluation, and averages the error over those samples.

Run the default study:

```bash
julia --project=. benchmarks/studies/gram_trajectory_lookahead_fidelity/main.jl
```

Useful options:

```bash
julia --project=. benchmarks/studies/gram_trajectory_lookahead_fidelity/main.jl \
  --planets mars,earth,venus,titan \
  --points 3,4,5,6,8,12,16,24,32,48,64 \
  --methods linear,pchip,cubic,akima \
  --horizon-s 900
```

Supported interpolation methods are `linear`, `pchip`, `cubic`, and `akima`.
`pchip`, `cubic`, and `akima` are implemented as cubic Hermite interpolants
with different slope estimates. The PCHIP-style option is shape-preserving and
usually less prone to overshoot than the unconstrained cubic option.

Outputs are written to `output/gram_trajectory_lookahead_fidelity/` by default:

- `gram_trajectory_lookahead_fidelity.csv`
- `gram_trajectory_lookahead_fidelity.pdf`
- `gram_trajectory_lookahead_fidelity.png`
- `gram_trajectory_lookahead_fidelity_relative.pdf`
- `gram_trajectory_lookahead_fidelity_relative.png`

The generated path is a simple linear altitude/latitude/longitude/elapsed-time
trajectory. Planet-specific default altitude ranges keep the query inside the
nominal atmospheric regions:

- earth: 130 km to 90 km
- mars: 145 km to 90 km
- venus: 170 km to 110 km
- titan: 700 km to 350 km

Override the shared path shape with `--h0-km`, `--h1-km`, `--lat0-deg`,
`--lat1-deg`, `--lon0-deg`, and `--lon1-deg`.

## Wind Interpolation Stress Test

Run a separate wind-focused stress test with broader altitude, latitude,
longitude, and elapsed-time trajectories:

```bash
julia --project=. benchmarks/studies/gram_trajectory_lookahead_fidelity/wind_interpolation_stress.jl
```

This writes `output/gram_wind_interpolation_stress/` by default:

- `gram_wind_interpolation_stress.csv`
- `gram_wind_interpolation_stress_rankings.csv`
- `gram_wind_interpolation_stress.pdf`
- `gram_wind_interpolation_stress.png`

The stress test ranks methods by wind-vector p95 error for each
planet/scenario/point-count combination.
