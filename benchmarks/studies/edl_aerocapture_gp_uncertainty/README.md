# EDL Aerocapture GP Uncertainty Study

This study evaluates whether a preceding aerocapture-like sensing pass can
improve real-time atmospheric uncertainty modeling for a later Earth EDL pass
relative to using Earth GRAM alone. The setup uses:

- a position-indexed truth field (see [Truth Sources](#truth-sources))
- a seeded `GRAMSuite.jl` nominal and ensemble as the onboard prior
- a Gaussian-process residual model over `lat, lon, alt` with an optional
  parametric mean function
- synthetic aerocapture and EDL corridors in nearby latitude/longitude bands

## Truth Sources

The estimator can only learn a residual that is a function of position. GRAM's
`perturbedDensity` is not: `PerturbedAtmosphere::setPosition` measures the
displacement since the previous call *on that instance* and accumulates it, so
the value returned depends on the instance's call history. Measured on this
study's own corridors:

- walking one instance along the same path twice: densities differ by up to `44.5%`
- aerocapture and EDL residuals at matched altitude in the `62-123 km` overlap:
  correlation `-0.13`, and `rms(edl - aero) = 0.1383` against `rms(edl) = 0.1365`,
  so subtracting the aerocapture measurement is worse than doing nothing

GRAM *nominal* density is deterministic in `(lat, lon, alt, t)` (repeatability
measured at exactly `0`), so every sound truth source is built on it.

| `--truth-source` | Truth | Use |
|---|---|---|
| `merra2_native` | Day-specific MERRA-2 **model-level** reanalysis, `72` levels to `0.01 hPa` (~`80 km`), `3-hourly` | The right answer: real weather, real time evolution, seam above the scored band. Needs a download |
| `merra2` (default) | The MERRA-2 grids vendored with Earth-GRAM, applied as a multiplicative anomaly on GRAM nominal | Climatology, ceiling ~`64 km`. Works with no download |
| `synthetic_field` | GRAM nominal times `exp(bias + amplitude * f)`, `f` an exact squared-exponential GP draw via random Fourier features | Validation: the planted signal is known, so a correct estimator must recover it |
| `gram_epoch_shift` | GRAM nominal at an epoch shifted from the one the prior assumes | Realistic same-model case: genuine day-to-day variability the prior cannot know |
| `gram_walk` | GRAM `perturbedDensity` walked along the trajectory | Legacy only. Not position-indexed; the driver warns and `RESULTS.md` is annotated. Retained to reproduce pre-fix numbers |

### MERRA-2

[`merra2.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/merra2.jl)
reads the binary grids vendored with Earth-GRAM directly, following the layout
in `Earth/source/MERRA2Data.cpp`: nineteen 3-D `(pressure, lat, lon)` Float64
blocks with density first and level height second, then twenty surface blocks.
The reader is validated three ways — the block count reproduces the exact
`106,012,800`-byte file size, `p / (rho T)` recovers `287.0 J/(kg K)` at every
probed level, and the interpolated profile matches GRAM's own nominal to within
`0.2%` between `20 km` and `45 km`, which is expected because GRAM builds its
nominal from the same data.

Truth is applied as a multiplicative anomaly so it stays defined at every
altitude:

```text
rho_truth = rho_gram * exp(taper * (log(rho_m2 / rho_gram) + dispersion * sigma_m2 * g))
```

Two properties of the vendored data drive the design:

- The grids are **monthly climatologies by time-of-day slot**, not specific-day
  reanalysis. The 18Z-minus-all-hours anomaly is only `0.05-0.5%` in
  log-density, an order of magnitude under the `5%` measurement noise, so a slot
  difference alone carries no usable signal. `--merra2-dispersion` therefore
  scales a smooth random field by MERRA-2's own relative interannual density
  standard deviation (`0.7-1.6%` on this corridor), giving a
  reanalysis-consistent specific day. Set it to `0` for the raw field.
- MERRA-2 spans `1000 mb` to `0.1 mb`, about `0.1 km` to **`64 km`**, not
  `80 km`. Above the local ceiling the anomaly is held at its ceiling value and
  tapered to zero over `--merra2-blend-width-km`, so truth reverts continuously
  to GRAM nominal — and the residual above roughly `63 km` is zero by
  construction.

That ceiling interacts directly with the corridor. See
[Sensing Depth](#sensing-depth).

## Mean Function

A zero-mean stationary GP decays to "no correction" more than a length scale
from its data. On the default corridor the maximum kernel correlation from an
EDL point near `35 km` to any aerocapture training point is `0.019`, so the
posterior mean there is `98%` prior by construction — the estimator cannot
differ from the prior in the band the metric weights most.

The dominant term turned out to be the assumed correlation length rather than
the geometry. See [Correlation Structure](#correlation-structure).

`--mean-basis` adds a parametric mean fitted by generalized least squares
(universal kriging), which extrapolates where the zero-mean GP term cannot, and
inflates the predictive variance where the basis itself is extrapolated:

- `none` — the original zero-mean residual GP
- `constant` — a bulk log-density offset relative to the prior
- `linear_alt` — offset plus an altitude trend

The fitted coefficients are written to `mean_function_fits.csv` and summarized
in `RESULTS.md`.

Which basis wins is diagnostic of the truth field's structure rather than a
fixed property of the estimator. On `synthetic_field` with a planted bulk offset
`constant` wins; on `merra2`, whose anomaly is spatially structured with a real
vertical profile, `none` wins at every sensing depth and `linear_alt` is badly
worse because it extrapolates a trend fitted from a thin sensed sliver.

## Correlation Structure

Earth-GRAM carries an explicit spatial correlation model. `EarthAtmosphere::getCorrelationCoefficients`
forms

```text
rho = exp(-( dx / Lh + dz / Lz + dt / Lt ))
```

a separable exponential on an L1 metric, with `Lh` and `Lz` tabulated against
altitude in `EarthAtmosphereData.cpp::rsData` and transcribed into
[`gram_correlation.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gram_correlation.jl):

| altitude | `Lh` | `Lz` |
|---|---|---|
| `20 km` | `2040 km` | `16.0 km` |
| `35 km` | `2279 km` | `18.3 km` |
| `50 km` | `2300 km` | `19.3 km` |
| `120 km` | `3900 km` | `32.0 km` |

### Is it present in MERRA-2?

Yes, horizontally, and to within a few percent through the band that matters.
[`measure_merra2_correlation.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/measure_merra2_correlation.jl)
builds two anomaly ensembles from the vendored grids over `0-45N` and fits an
`e`-folding scale to the binned correlation against chordal separation:

| altitude | MERRA-2 `e`-fold | GRAM `Lh` | ratio |
|---|---|---|---|
| `8.4 km` | `6048 km` | `821 km` | `7.37` |
| `16.4 km` | `1657 km` | `1582 km` | `1.05` |
| `26.3 km` | `1263 km` | `2156 km` | `0.59` |
| `35.6 km` | `2663 km` | `2283 km` | `1.17` |
| `39.3 km` | `2333 km` | `2305 km` | `1.01` |
| `50.4 km` | `2682 km` | `2295 km` | `1.17` |
| `56.7 km` | `1353 km` | `2214 km` | `0.61` |

Two caveats on the method. The horizontal ensemble is the eight time-of-day
slots per month as departures from that month's all-hours mean, so it measures
the correlation structure of the *diurnal* signal — broad and tide-driven, which
is why the `8.4 km` row is an outlier rather than a refutation of GRAM's `821 km`
there. And the vertical structure cannot be read from that ensemble at all,
because migrating tides propagate vertically with a wavelength around `20-30 km`
and make the correlation oscillate in sign instead of decaying. The vertical
column therefore uses the twelve monthly means as departures from the annual
mean, which is deep and coherent but only twelve members. On that ensemble the
decay about `30.8 km` tracks GRAM's `exp(-|dz| / 17.4)` within roughly `0.1` out
to `+16 km` and is noisy below.

### Does it improve extrapolation?

Yes, but the effect is second-order once a mean basis is fitted, and the tables
that used to sit here were measured against the superseded climatology truth.
See [Current Results](#current-results) for the numbers that stand.

What survives from that work is the *scale* measurement, which is independent of
which truth source is used: the study's original `1.5 deg` (about `166 km`)
horizontal length scale was roughly `14x` shorter than both GRAM's tabulated
`Lh` and MERRA-2's measured structure. The CLI defaults are now GRAM's values
(`--lat-scale-deg 20.5 --lon-scale-deg 21.8 --alt-scale-km 18.33`), and
`gram_exponential` uses GRAM's own separable-exponential form directly.

### Time

GRAM's correlation includes a temporal term, `dt / Lt` with
`Lt = max(3 hr, 0.735 day * h_km^0.116)` — about `26.6 hr` at `35 km`, so a
`6 hr` gap retains `rho_t = 0.80`. `--kernel-time-scale-hours` adds the matching
factor to the kernel.

Measured against real reanalysis, the gap costs far less than that model implies:
`45.4% / 44.3% / 45.1%` at `1 / 3 / 6 hr` on the `75-35 km` band, essentially
flat. `Lt` is GRAM's *small-scale perturbation* timescale, and a synoptic-scale
density anomaly persists considerably better over six hours than it suggests.

The `--merra2-time-decorrelation` AR(1) surrogate exists only for the
climatology truth source, which contains no weather and therefore no real
evolution. `merra2_native` supplies genuine `3-hourly` evolution and does not
use it.

## Uncertainty: `--amplitude-mode`

The original estimator gave the GP one global kernel amplitude,
`var(detrended residuals)`, and then re-introduced the GRAM prior variance
afterwards with

```text
fused_mu  = prior_mu + delta_mu                        # additive
fused_var = 1/(1/prior_var + 1/delta_var)              # inverse-variance
```

Those come from two different models. The additive form pairs with
`var = prior_var + delta_var`; the inverse-variance form pairs with a
precision-weighted mean. Taking one from each means the mean moves by the full
GP correction while the variance is computed as if the GP were a redundant
second measurement of the same quantity. Since `1/prior_var > 0`, the result is
always smaller than `delta_var`, and for a zero-mean GP `delta_var <= amplitude^2`
— so the reported variance could never exceed the kernel amplitude anywhere in
the corridor, however uncertain the prior or however far from data the point.

`--amplitude-mode prior_scaled`, the default, folds the prior variance into the
kernel instead:

```text
k(x, x') = lambda^2 * s(x) * s(x') * rho(x, x')
```

`s(x)` is the GRAM ensemble relative sigma and `rho` the correlation function;
`D K D` with `D = diag(s)` is positive definite whenever `K` is. The GP prior
variance at `x` is then exactly `lambda^2 s(x)^2`, so the posterior returns to
the GRAM prior far from data and tightens only near it, and the predictive
variance *is* the GP posterior variance — there is no fusion step to get wrong.
`lambda` is a single inflation factor for GRAM ensemble mis-dispersion, fitted
by maximizing the restricted log marginal likelihood.

This changes the mean as well. In the noise-free limit,

```text
mu(x*) = s(x*) * interp( y_i / s(x_i) )
```

so residuals transfer in units of the local prior sigma rather than in absolute
log-density. "The atmosphere is `+1.5 sigma`" propagates, instead of "`+13%` at
`90 km`" propagating unchanged to `35 km` where the real spread is `1.2%`.
`--amplitude-mode stationary` reproduces the legacy behaviour.

## Scope

The implemented study makes the favorable assumptions requested for a
preliminary test:

- aerocapture and EDL stay near the same ground track
- the time gap is a few hours
- trajectory knowledge is perfect
- aerocapture yields direct noisy density samples every `2 s`
- measurement noise is `N(0, (0.05 rho_truth)^2)`
- the GP is held fixed after the aerocapture pass
- there is no onboard density filter during EDL; the baseline is GRAM alone
- the truth anomaly decorrelates over the gap on GRAM's own timescale, and the
  kernel knows it (see [Time](#time))

The estimator compares three state parameterizations:

- `log_density`
- `density_scale_factor`
- `log_density_scale_factor`

and four kernels:

- squared exponential
- Matérn `3/2`
- Matérn `5/2`
- `gram_exponential`, GRAM's own separable-exponential correlation model on
  physical distance with its tabulated altitude-dependent scales

## Layout

- [`main.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/main.jl): study driver
- [`corridor.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/corridor.jl): synthetic aerocapture and EDL paths
- [`truth_sources.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/truth_sources.jl): pluggable truth fields
- [`merra2.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/merra2.jl): reader for the vendored MERRA-2 binary climatology grids
- [`merra2_native.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/merra2_native.jl): reader for day-specific MERRA-2 model-level netCDF granules
- [`fetch_merra2_native.sh`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/fetch_merra2_native.sh) and [`verify_merra2_native.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/verify_merra2_native.jl): download and post-download checks
- [`gram_correlation.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gram_correlation.jl): GRAM's tabulated correlation scales and the chordal metric
- [`prior_sources.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/prior_sources.jl): swappable onboard prior (GRAM nominal or NRLMSISE-00)
- [`accel_filter.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/accel_filter.jl) and [`compare_filter.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/compare_filter.jl): in-situ density-ratio filter and its comparison against the GP
- [`measure_merra2_correlation.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/measure_merra2_correlation.jl): measures MERRA-2's empirical correlation structure
- [`gram_prior.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gram_prior.jl): seeded GRAM nominal and dispersed prior sampling
- [`gp_models.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gp_models.jl): GP residual models and target transforms
- [`scoring.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/scoring.jl): RMSE, NLPD, and coverage metrics
- [`plot_results.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/plot_results.jl): summary plotting
- [`RESULTS.md`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/RESULTS.md): generated summary after a run

## Project Setup

This study uses its own Julia environment for `SpaceAGORA.jl` and
`GRAMSuite.jl`.

Instantiate it with:

```text
julia --project=benchmarks/studies/edl_aerocapture_gp_uncertainty -e "using Pkg; Pkg.instantiate()"
```

The study project declares:

- `SpaceAGORA` from the local repository root
- `GRAMSuite` from `data/GRAMSuite.jl`

You still need:

- a working local GRAM asset install under `data/GRAMSuite.jl/GRAM Suite 2.0`
- a usable native GRAM build
- the MERRA-2 grids under `.../Earth/data/MERRA2data`, or `SPACEAGORA_MERRA2_PATH`
  pointing at them, for the default `merra2` truth source

## Running

Default run:

```text
julia --project=benchmarks/studies/edl_aerocapture_gp_uncertainty benchmarks/studies/edl_aerocapture_gp_uncertainty/main.jl
```

Useful options:

```text
--out-dir output/edl_aerocapture_gp_uncertainty
--seed 42
--truth-seed 10042
--prior-seed 20042
--n-dispersion 24
--case-limit 1
--save-pointwise true
--lat-scale-deg 20.5
--lon-scale-deg 21.8
--alt-scale-km 18.33
--mean-basis none,constant,linear_alt
--kernels squared_exponential,matern32,matern52,gram_exponential
--gram-scale-ref-alt-km 35.0
--kernel-time-scale-hours auto
--amplitude-mode prior_scaled
--merra2-time-decorrelation true
--aerocapture-exit-alt-km 60.0
--truth-source merra2
--merra2-native-dir <path>
--merra2-hour-code 0
--merra2-dispersion 1.0
--merra2-blend-width-km 10.0
--truth-epoch-shift-hours 36.0
--field-bias 0.08
--field-amplitude 0.10
--field-lat-scale-deg 6.0
--field-lon-scale-deg 12.0
--field-alt-scale-km 25.0
```

## Outputs

The driver writes:

- `summary_metrics.csv`
- `kernel_comparison.csv`, averaged by parameterization, kernel, and mean basis
- `mean_function_fits.csv`, the fitted GLS mean coefficients per case
- `cases.csv`, including the truth source and whether it is position-indexed
- `pointwise_predictions.csv` when enabled
- a generated [`RESULTS.md`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/RESULTS.md)

`plot_results.jl` reads the summary and pointwise CSV files and produces bar
charts for weighted metrics plus a representative density-profile comparison
and relative-error profile. Pointwise output must be enabled for the profile
plots. Generated figures are saved under `plots/` in this study directory by
default.

## Method Notes

The onboard prior uses the nominal GRAM density with `--prior-seed`; its
dispersed ensemble uses GRAM's seeded `perturbedDensity` with consecutive,
different seeds. The defaults derived from `--seed` are distinct
(`seed + 10_000` and `seed + 20_000`) and are recorded in `cases.csv`. This is
a synthetic same-model experiment, not an independent atmospheric-truth
validation.

The GP is implemented as a residual model on top of the GRAM prior mean. The
GRAM dispersed ensemble is used to construct the baseline uncertainty and a
pointwise prior variance for the transformed variables. The fusion of prior and
GP variances is intentionally simple and should be treated as a first-pass
uncertainty fusion heuristic rather than a final flight estimator design.

The EDL weighting profile emphasizes altitudes near a nominal Earth peak
dynamic-pressure region, implemented here as a Gaussian weight centered near
`35 km`.

## `merra2_native`: Day-Specific Model-Level Reanalysis

The seam above is a coverage problem. GRAM ships the MERRA-2 **pressure-level**
grid — Table 4.1 of the MERRA-2 File Specification, `1000 hPa` to `0.1 hPa`,
exactly the `42` levels in `merra2.jl` — averaged into monthly time-of-day
climatologies. Two consequences followed from that and both are fixed by
switching products rather than by any modelling change:

| | `merra2` (vendored climatology) | `merra2_native` (M2I3NVASM) |
|---|---|---|
| vertical grid | `42` pressure levels | `72` hybrid sigma-p model levels |
| ceiling | `0.1 hPa`, about `64 km` | `0.01 hPa` (PTOP), about `80 km` |
| time | monthly climatology by time-of-day slot | `3-hourly` analysis, per day |
| weather | none — averaged out | real |
| level order | bottom-up | **top-down**, level 1 is the model top |

With `--truth-source merra2_native` the synthetic dispersion term and the AR(1)
time surrogate both become unnecessary: the anomaly is real weather and it
evolves genuinely across the aerocapture-to-EDL gap. The ceiling moves about
`16 km`, putting the GRAM/MERRA-2 blending seam well above the `24-46 km`
scored band.

### Getting the data

Granules are not vendored — `~2.1 GB` each, and GES DISC requires a NASA
Earthdata login:

```text
./benchmarks/studies/edl_aerocapture_gp_uncertainty/fetch_merra2_native.sh \
    2024-03-20 2024-03-21 2024-10-05 2024-10-06 2025-03-21 2025-03-22
```

A `6 hr` gap from an `18:00 UTC` anchor crosses midnight, which is why each
anchor needs the following day too. `--dry-run` prints the URLs without
fetching, and the script explains the Earthdata authorization step that causes
most `401`s. GES DISC's spatial subsetter cuts a corridor-sized window to a few
tens of MB per day if bandwidth is tight.

Then, before trusting any result:

```text
julia --project=benchmarks/studies/edl_aerocapture_gp_uncertainty \
      benchmarks/studies/edl_aerocapture_gp_uncertainty/verify_merra2_native.jl
```

### What is verified, and what is not

The reader's logic is pinned by `runtests.jl`, which writes a synthetic granule
carrying the documented schema — `tzyx` variables presented as
`(lon, lat, lev, time)`, `72` top-down levels, `_FillValue = 1e15` — with a
planted atmosphere whose log-density is linear in latitude, longitude, time and
altitude. Every interpolation the reader performs is exact on that field, so
recovery is checked to `1e-6` relative. That pins the level flip, the fill
handling, the bilinear and log-linear interpolations and the time axis.

Two things cannot be settled without real granules, and
`verify_merra2_native.jl` reports on both:

- **Whether `H` is geopotential or geometric height.** The reader assumes
  geopotential and converts, a `~1 km` difference at `80 km` — not negligible
  against an `18 km` vertical correlation scale. Check 3 prints both readings
  against GRAM's geometric-altitude profile; the wrong one shows a systematic
  ramp with altitude.
- **Whether a subsetted file kept its level ordering.** The specification warns
  that some post-processors flip the vertical grid. The reader asserts the level
  count, and check 2 confirms density falls with altitude.

## Current Results

All numbers below are from `output/edl_gp_final` and `output/edl_gp_e75`: 36
cases (3 anchors x 3 gaps x 4 ground-track offsets), truth `merra2_native`,
prior `gram`, `--amplitude-mode prior_scaled`. Regenerate with
`--report-only true`; `RESULTS.md` and `plots/` track the most recent run.

Nominal coverage is `0.683` and `0.954`. The GRAM baseline scores
`0.03398` weighted log RMSE at `0.350` / `0.405` coverage.

| sensing band | mean basis | vs GRAM | weighted NLPD | cov `1-sigma` | cov `2-sigma` |
|---|---|---|---|---|---|
| `125-35 km` | `constant` | `39.53%` | `-5.01` | `0.364` | `0.555` |
| `125-35 km` | `linear_alt` | `44.74%` | `-8.79` | `0.456` | `0.725` |
| `125-35 km` | `none` | `23.26%` | `-6.69` | `0.366` | `0.484` |
| **`75-35 km`** | **`constant`** | **`47.93%`** | **`-8.57`** | **`0.559`** | **`0.768`** |
| `75-35 km` | `linear_alt` | `22.58%` | `-8.53` | `0.464` | `0.771` |
| `75-35 km` | `none` | `37.25%` | `-8.47` | `0.507` | `0.738` |

`gram_exponential` with a `constant` mean over a `75-35 km` sensing band is the
best configuration on both axes. Kernel choice is second-order once a mean basis
is fitted; all four land within a few points of each other.

### Why the sensing band matters more than any estimator choice

Lowering the entry from `125 km` to `75 km` is worth more than any kernel or
parameterization change. Above the `52 km` anomaly cap the truth residual is
zero by construction, so a pass starting at `125 km` spends most of its length
measuring nothing but the `5%` sensor noise. That is directly visible in the
fitted inflation factor:

| sensing band | share of pass with signal | fitted `lambda` (`none`) |
|---|---|---|
| `125-35 km` | `18.9%` | `1.00`, pinned at the floor |
| `75-35 km` | `42.5%` | `1.52`, freely identified |
| `55-35 km` | `85.0%` | `3.02` |

`linear_alt` scoring `44.74%` at `125 km` and `22.58%` at `75 km` is not a
property of the estimator: over the longer pass the linear trend is anchored by
a large block of zero-residual points above `52 km`, and it was partly fitting
that artifact.

`signal_fraction` stays at `0.09-0.14` throughout. Only about an eighth of the
training-residual variance sits above the measurement-noise floor, and every
result here rests on that.

### Comparison with an in-situ density filter

[`compare_filter.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/compare_filter.jl)
runs a scalar density-ratio filter in the style of Tracy, Falcone and
Manchester, *Robust Entry Guidance with Atmospheric Adaptation*, against the
same truth. Twelve cases, MSL-class vehicle (`beta = 121 kg/m^2`), corridor
speed `4.0 km/s`:

| | at the vehicle | vs prior | forward | vs prior |
|---|---|---|---|---|
| prior alone | `0.03814` | — | `0.04148` | — |
| in-situ filter only | `0.02305` | `39.5%` | `0.04934` | **`-19.0%`** |
| aerocapture GP only | `0.01959` | `48.6%` | **`0.01679`** | **`59.5%`** |
| GP + filter, naive product | `0.02363` | `38.1%` | `0.04212` | `-1.5%` |
| filter, correlation-decayed | `0.02305` | `39.5%` | `0.03893` | `6.2%` |
| GP + filter, correlation-decayed | `0.02363` | `38.1%` | `0.01941` | `53.2%` |

"forward" is the error over the *remaining* trajectory as predicted at each
guidance call, which is what a predictor-corrector consumes: CPEG re-simulates
to parachute deploy every call.

The filter is strong where it is measuring and **worse than doing nothing ahead
of the vehicle**: a scalar held forward smears a local ratio across a region
where the true ratio varies with altitude. Decaying that correction over GRAM's
vertical correlation scale recovers it to `6.2%`, still far below the GP's
`59.5%`.

Naively multiplying the two is worse than either alone, for the same reason.
Decaying the coupling fixes that (`53.2%`), but it does not beat the GP by
itself here — the GP already carries the spatial structure the filter lacks, and
a scalar correction on top mostly adds noise.

Two caveats on this comparison. The filter is driven by **drag acceleration
directly**, where the paper uses navigation position and velocity; a direct
acceleration measurement is strictly more informative about instantaneous
density, so these filter numbers are an upper bound on that architecture. And
this study's corridor is kinematic — it does not respond to density — so nothing
here is closed-loop. The paper's actual result is a `37%` reduction in mean
parachute-deployment miss distance, which an open-loop density metric cannot
speak to.

## Known Remaining Issues

These are diagnosed but not yet fixed, and they bound what the current numbers
can show:

- The aerocapture pass bottoms out at `60 km` while the `q`-weight profile is
  centered at `35 km`, so the sensed and scored bands do not overlap.
  `--aerocapture-exit-alt-km` now exposes this and
  [Sensing Depth](#sensing-depth) quantifies it, but the default corridor is
  still mismatched.
- The MERRA-2 truth residual is zero above roughly `63 km` by construction,
  because MERRA-2 has no data there. The upper part of the sensing pass
  therefore carries no information under that truth source. An upper-atmosphere
  truth field would need a different data source.
- `gram_exponential` freezes GRAM's altitude-dependent scales at a single
  reference altitude, because a kernel with position-dependent scales is not
  positive definite without a Gibbs-style normalization that has no clean form
  for the exponential. Over `30-60 km` GRAM's `Lh` spans `2240-2318 km` and `Lz`
  spans `17.4-19.3 km`, so the approximation is mild, but it is an approximation.
- The AR(1) time decorrelation applies to the *whole* anomaly, including the
  `log(rho_m2 / rho_gram)` model-difference component. A genuine systematic bias
  between GRAM and MERRA-2 would persist rather than decay, so with
  `--merra2-dispersion 0` the construction is conservative: it decorrelates
  something that in reality would partly hold. Splitting persistent model bias
  from decaying weather would need a truth source that separates them.
- `Lt` is GRAM's scale for its *small-scale* perturbation model. A synoptic
  density anomaly plausibly persists longer and a gravity-wave component
  decorrelates in well under an hour, so a single `Lt` is a coarse stand-in for
  a spectrum of timescales.
- Kernel hyperparameters are set from the CLI and the amplitude from the
  detrended sample variance; nothing is fitted by marginal likelihood.
- The `52 km` anomaly cap means a pass entering at `125 km` spends most of its
  length measuring noise. That is a truth-construction artifact, not physics,
  and it is what drives `lambda` to the floor. Filling `52-125 km` with a real
  truth source (NAVGEM-HA, SD-WACCM-X or SABER) would remove it; lowering the
  entry to `75 km` sidesteps it at the cost of scenario realism.
- Nothing here is closed-loop. The comparison against the Tracy/Falcone/
  Manchester filter scores open-loop density prediction; their reported result
  is a `37%` reduction in parachute-deployment miss distance, which this
  framework cannot measure.
- `lambda` is fitted on the sensed band and applied to the scored band. On this
  corridor the two disagree: the implied inflation from the residuals is about
  `1.6` over `45-63 km` and about `0.47` over `25-45 km`. A single scalar cannot
  reconcile them, and no aerocapture pass can measure GRAM's ensemble
  calibration at an altitude it never reaches.
- `weighted_rmse` is an absolute-density metric, and on the default corridor the
  `10-20 km` band carries `87%` of its weighted squared error while the `q`-max
  band carries `0.8%`. Use `weighted_rmse_log`. The `q` weights themselves are a
  Gaussian in altitude with a `0.25` floor rather than a real `0.5 rho V^2`
  profile.
- `log_density` and `log_density_scale_factor` are the same estimator and score
  identically to `14` digits.

## Verification Status

The seed-plan and driver-load smoke tests pass. The first native-GRAM
validation step is a reduced single-case run with
`--n-dispersion 4` and `--save-pointwise false` to confirm:

- the truth, nominal, and ensemble GRAM seeds produce distinct realizations
- GRAM nominal and dispersed point sampling work
- the generated `summary_metrics.csv` looks numerically sane
