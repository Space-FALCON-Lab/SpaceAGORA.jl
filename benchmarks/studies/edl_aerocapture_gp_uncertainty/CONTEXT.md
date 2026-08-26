# Conversation Context

## Study Goal

This study is intended to test whether an aerocapture-like sensing pass can
improve atmospheric uncertainty modeling for a later EDL pass relative to using
GRAM alone.

The requested preliminary-test assumptions are:

- Earth is used as a proxy environment.
- A separately seeded Earth GRAM realization provides temporary synthetic truth.
- Earth GRAM nominal and a separately seeded dispersed ensemble provide the onboard prior.
- Aerocapture and EDL occur in approximately the same latitude/longitude band.
- The time gap between aerocapture and EDL is a few hours.
- The GP should infer over `lat, lon, alt`.
- The time component was deferred for now; if there are no new measurements,
  the GP is held fixed.
- Synthetic aerocapture density samples are taken every `2 s`.
- Measurement noise is Gaussian with zero mean and standard deviation equal to
  `5%` of truth density.
- Perfect trajectory knowledge is assumed.
- Aerocapture starts at `125 km` and descends to about `60 km`.
- EDL starts at `125 km` and descends to about `10 km`.
- Performance is assessed on atmospheric prediction quality only, not closed-loop guidance.
- Metrics should emphasize the region near nominal peak dynamic pressure.

## Requested Comparisons

The study should compare:

- `log_density`
- `density_scale_factor`
- `log_density_scale_factor`

and these kernels:

- squared exponential
- Matérn `3/2`
- Matérn `5/2`

The GRAM prior contribution should use:

- GRAM nominal as the prior mean
- a variance estimated from dispersed GRAM samples

The baseline to beat is GRAM predictive performance including uncertainty, not
just the GRAM mean.

## Implemented Files

The following scaffold was added under this study directory:

- [`Project.toml`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/Project.toml)
- [`main.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/main.jl)
- [`corridor.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/corridor.jl)
- [`truth_sources.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/truth_sources.jl)
- [`merra2.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/merra2.jl)
- [`gram_correlation.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gram_correlation.jl)
- [`measure_merra2_correlation.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/measure_merra2_correlation.jl)
- [`truth_cache.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/truth_cache.jl)
- [`gram_prior.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gram_prior.jl)
- [`gp_models.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/gp_models.jl)
- [`scoring.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/scoring.jl)
- [`plot_results.jl`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/plot_results.jl)
- [`README.md`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/README.md)
- [`RESULTS.md`](/home/space-falcon-1/Documents/SpaceAGORA.jl/benchmarks/studies/edl_aerocapture_gp_uncertainty/RESULTS.md)

## Diagnosis of the Weak First Result

The first full run improved weighted RMSE over the GRAM baseline by `0.9%` and
scored worse `2-sigma` coverage than the baseline. Three independent causes were
identified and measured:

1. GRAM `perturbedDensity` is a call-history random walk, not a spatial field.
   One instance walked twice along the same path differs by up to `44.5%`, and
   the aerocapture and EDL residuals at matched altitude correlate at `-0.13`.
   No position-indexed estimator can learn it.
2. The GP was extrapolating `2-4` length scales. Max kernel correlation from an
   EDL point at `33 km` to any training point is `0.019`, so the posterior mean
   was `98%` prior by construction.
3. `weighted_rmse` is absolute-density RMSE, and the `10-20 km` band carries
   `87%` of its weighted squared error while the `q`-max band carries `0.8%`.

Two bugs were found alongside these: the precision-weighted variance fusion
shrinks uncertainty where the GP has no information, and `log_density` and
`log_density_scale_factor` are the same estimator.

## Fixes Applied

Truth sources are now pluggable through `truth_sources.jl`, with
`synthetic_field` (validation, planted signal) and `gram_epoch_shift`
(realistic) replacing the unsound `gram_walk` as defaults. The residual GP
gained an optional GLS mean function (`--mean-basis none|constant|linear_alt`),
which extrapolates a bulk density offset into the unsensed band, and
`weighted_rmse_log` was added so an improvement in the sensed band is visible.

Measured on `synthetic_field` with a planted `bias = 0.08`, three cases,
`--n-dispersion 8`:

| mean basis | weighted log RMSE vs GRAM | fitted `beta_constant` |
|---|---|---|
| `none` | `37.2%` | n/a |
| `constant` | `52.6%` | `0.0785-0.0810` |
| `linear_alt` | `-3.5%` | `0.0789-0.0810` |

On `gram_epoch_shift` the zero-mean GP reaches `64.5%` and the constant basis
`39.6%`: the realistic residual is a smooth vertical structure rather than a
bulk offset, so forcing a constant biases the extrapolation. `linear_alt` is
worse than `constant` on the synthetic field because that field has no genuine
altitude trend, so the linear term extrapolates noise.

## MERRA-2 Truth Source

`merra2.jl` reads the binary grids vendored with Earth-GRAM directly (layout
from `Earth/source/MERRA2Data.cpp`: nineteen 3-D Float64 blocks, density first,
level height second, then twenty surface blocks). Validated by exact file-size
arithmetic, by `p / (rho T) = 287.0`, and by matching GRAM's own nominal to
`0.2%` between `20` and `45 km`.

Two data properties shaped the design. The grids are monthly climatologies by
time-of-day slot, not specific-day reanalysis: the 18Z-minus-all-hours anomaly
is only `0.05-0.5%` in log-density, well under the `5%` measurement noise, so
`--merra2-dispersion` scales a smooth field by MERRA-2's own relative
interannual sigma (`0.7-1.6%` here) to build a reanalysis-consistent day.
And MERRA-2 tops out at `0.1 mb`, about `64 km`, not the `80 km` assumed
earlier, so the anomaly is tapered back to GRAM nominal above the local ceiling
and the residual above roughly `63 km` is zero by construction.

That ceiling turned out to dominate the result. Sweeping
`--aerocapture-exit-alt-km`, three cases, `--n-dispersion 8`, weighted log RMSE
against a `0.030393` baseline:

| exit altitude | `none` | `constant` | `linear_alt` |
|---|---|---|---|
| `60 km` | `2.95%` | `-1.92%` | `-202.9%` |
| `55 km` | `3.65%` | `-0.97%` | `-115.0%` |
| `45 km` | `3.21%` | `-5.33%` | `-61.2%` |
| `35 km` | `31.50%` | `25.19%` | `2.32%` |
| `25 km` | `54.50%` | `48.99%` | `43.13%` |

The jump between `45` and `35 km` tracks whether the sensed band reaches the
scored band (the `q` weight spans about `24-46 km`), not how many training
points fall inside MERRA-2's domain, which grows smoothly from `23%` to `32%`
across that step. A realistic Earth aerocapture does not reach `25 km`, so
either the pass must go deeper than an aerocapture plausibly would or the study
must score the band it actually senses.

On MERRA-2 truth `none` beats both mean bases at every depth: the anomaly is
spatially structured with a real vertical profile rather than a bulk offset, so
the parametric mean fits a thin sensed sliver and extrapolates it badly. Which
basis wins is diagnostic of the truth field, not a fixed property of the
estimator.

## Spatial Correlation

Earth-GRAM's `EarthAtmosphere::getCorrelationCoefficients` uses
`rho = exp(-(dx/Lh + dz/Lz + dt/Lt))`, a separable exponential on an L1 metric,
with `Lh`/`Lz` tabulated in `EarthAtmosphereData.cpp::rsData`
(`Lh = 2279 km`, `Lz = 18.3 km` at `35 km`). Transcribed into
`gram_correlation.jl` and added as the `gram_exponential` kernel.

MERRA-2 confirms the horizontal structure to within a few percent between
`16 km` and `50 km` (`measure_merra2_correlation.jl`): measured `e`-folding
`1657/1263/2663/2333/2682 km` at `16.4/26.3/35.6/39.3/50.4 km` against GRAM's
`1582/2156/2283/2305/2295 km`. The `8.4 km` outlier (`7.4x`) is the
tide-dominated troposphere in the diurnal ensemble, not a refutation. The
vertical structure needs the seasonal ensemble instead, because migrating tides
propagate with a `20-30 km` wavelength and make the diurnal vertical correlation
oscillate in sign; on the seasonal ensemble the decay tracks GRAM's
`exp(-|dz|/17.4)` within about `0.1` out to `+16 km`.

This revises the earlier diagnosis. Cause 2 was attributed to corridor geometry,
but the dominant term is the assumed correlation length. Truth `merra2`, three
cases, mean basis `none`, vs a `0.030393` baseline:

| exit | old scales `matern32` | GRAM scales `matern32` | `gram_exponential` |
|---|---|---|---|
| `60 km` | `2.95%` | `-26.17%` | `12.76%` |
| `55 km` | `3.65%` | `44.96%` | `51.46%` |
| `45 km` | `3.21%` | `56.96%` | `57.14%` |
| `35 km` | `31.50%` | `53.92%` | `55.57%` |

Isolated at `45 km` exit with `matern32`: horizontal `166 km` + vertical `12 km`
gives `3.21%`; `166 km` + `18.33 km` gives `-8.61%`; `2276 km` + `12 km` gives
`64.61%`; `2276 km` + `18.33 km` gives `56.96%`. The horizontal scale is the
whole effect, and lengthening the vertical alone is actively harmful because a
short horizontal scale plus a long vertical one just propagates local noise down
the column. The `1.5 deg` default was `14x` too short; CLI defaults are now
GRAM's values.

The `60 km` row is instructive: with the pass almost entirely above MERRA-2's
ceiling there is nearly no signal, and a smooth long-scale `matern32` spreads a
fit of noise across the domain and loses to the prior, while the rougher
`gram_exponential` reverts faster and stays positive.

## Time Decorrelation

The study was measuring a frozen atmosphere: `truth_profiles(::Merra2Truth, ...)`
loaded one MERRA-2 grid and evaluated it with no time argument, so the anomaly
at aerocapture was bit-identical to the anomaly at EDL hours later. The baseline
error was the same to five significant figures across `1/3/6 hr` gaps.

Fixed with GRAM's own temporal scale, `Lt = max(3 hr, 0.735 day * h_km^0.116)`,
about `26.6 hr` at `35 km`. Truth evolves as a stationary AR(1),
`a(x,t) = rho_t a(x,0) + sqrt(1-rho_t^2) sigma_a g_innov(x)`, with `g_innov`
carrying GRAM's spatial scales; `--kernel-time-scale-hours` adds the matching
`dt/Lt` factor to the kernel. Both default on.

`--merra2-dispersion 0`, `45 km` exit, `gram_exponential`, twelve cases:

| truth | GP time axis | `1 hr` | `3 hr` | `6 hr` |
|---|---|---|---|---|
| frozen | off | `52.12%` | `40.26%` | `42.93%` |
| AR(1) | off | `45.00%` | `23.81%` | `13.48%` |
| AR(1) | GRAM `Lt` | `44.84%` | `26.50%` | `19.75%` |

The frozen field inflated the `6 hr` result about `3x` and produced no gap
dependence; the time-aware kernel recovers more of the loss the larger the gap.
Calibration at `6 hr` improves too, NLPD `-7.705` to `-8.193` and `2-sigma`
coverage `0.516` to `0.598`, but both stay well short of the baseline's `-8.967`
and `0.884` because the precision-weighted variance fusion is still unfixed.

Two caveats. The AR(1) decorrelates the whole anomaly including the
GRAM-vs-MERRA-2 model difference, which in reality would partly persist, so with
`dispersion 0` it is conservative. And `Lt` is GRAM's small-scale perturbation
timescale, a coarse stand-in for a spectrum from gravity waves to synoptic.

## Variance Fusion Replaced

`predict_density_profile` combined an additive mean with an inverse-variance
fusion — two different models. Because `1/prior_var > 0`, the result was always
below `delta_var`, and for a zero-mean GP `delta_var <= amplitude^2`, so the
reported variance could never exceed the kernel amplitude anywhere.

Replaced by `--amplitude-mode prior_scaled` (now the default):
`k(x,x') = lambda^2 s(x) s(x') rho(x,x')` with `s(x)` the GRAM ensemble relative
sigma, `lambda` fitted by restricted marginal likelihood. The GP prior variance
is then exactly `lambda^2 s(x)^2`, the posterior returns to the GRAM prior far
from data, and the predictive variance is the posterior variance with no fusion.
It also changes the mean: `mu(x*) = s(x*) interp(y_i / s(x_i))`, so residuals
transfer in units of the local prior sigma instead of absolute log-density.
`--amplitude-mode stationary` reproduces the legacy path.

`gram_exponential`, mean basis `none`, three cases; nominal coverage `0.683`/`0.954`:

| config | metric | stationary | prior_scaled | baseline |
|---|---|---|---|---|
| `60 km` exit, disp 1 | weighted log RMSE | `0.0300` | `0.0200` | `0.0299` |
| | weighted NLPD | `-7.293` | `-9.568` | `-9.502` |
| | weighted cov `1sig` | `0.275` | `0.744` | `0.880` |
| | weighted cov `2sig` | `0.535` | `0.928` | `0.921` |
| `45 km` exit, disp 0 | weighted log RMSE | `0.0197` | `0.0215` | `0.0334` |
| | weighted NLPD | `-9.292` | `-9.478` | `-9.299` |
| | weighted cov `1sig` | `0.420` | `0.687` | `0.801` |
| | weighted cov `2sig` | `0.855` | `0.926` | `0.900` |

Predictions made before implementing were wrong in both directions and the
reason is worth keeping. The fusion's shrink factor was measured as a near
no-op in the `25-45 km` band *in the default configuration*, and that was
generalized into "weighted calibration will not move" and "RMSE will move a few
points either way". Both failed: weighted `1-sigma` coverage went `0.275` to
`0.744`, weighted NLPD closed the gap to the baseline and passed it, and
weighted RMSE went from no improvement at all to `33%`. The error was treating a
variance-side measurement in one configuration as if it bounded the *mean*-side
effect in another. The mean change — transferring normalized rather than
absolute residuals — was the dominant term.

`lambda` fits to about `1.04`, not the `1.6` estimated from raw z-scores, because
marginal likelihood correctly attributes most of the sensed-band residual to the
`5%` measurement noise.

## Headline Numbers Are Seam-Dominated

Band decomposition of the default run's weighted log-density error: `55-73 km`
carries `94%` of the GRAM baseline error, and that band is where GRAM's nominal
blends away from MERRA-2 near the `0.1 mb` ceiling (raw
`log(rho_m2/rho_gram)` goes `-0.011` at `52.7 km` to `-0.170` at `62.6 km`) plus
this study's own taper ramp. In the `q`-max band `24-46 km` the GP is `-135%`
against GRAM alone at a `60 km` exit and `-117%` at `45 km`; it also adds error
above `73 km` where truth equals the prior by construction.

The GP helps in and just below the sensed band and degrades below it. All
whole-corridor improvement percentages reported earlier are largely seam
learning. Fix by tapering from well below the ceiling or scoring below `55 km`,
then re-run everything.

An oracle scalar density ratio, the quantity Tracy/Falcone/Manchester's SQEKF
estimates, explains only `6.9%` of the log-density error in the `q`-max band and
`1.4%` whole-corridor here — yet their closed-loop result is a `37%` mean
miss-distance reduction. Pointwise density RMSE is a far harsher scoreboard than
guidance actually needs.

## Current State

The study now uses independent GRAM seed streams for synthetic truth and the
GRAM prior. Julia is executable in this shell and the study project
instantiates successfully.

Remaining, in priority order:

1. Couple the GP to a conventional accelerometer-based density filter during
   EDL, and rebaseline against `GRAM + filter` rather than GRAM alone.
2. Move to closed-loop targeting-error scoring, since lead time does not show up
   in an open-loop pointwise density metric.
3. Decide the corridor question the sweep exposes: fly the sensing pass deeper,
   or score the band it actually senses.
4. Fit kernel hyperparameters by marginal likelihood, now that the defaults are
   at least physically motivated.
5. Reconcile `lambda` across bands, or accept that the aerocapture pass cannot
   calibrate GRAM's ensemble at altitudes it never reaches.
6. Find an upper-atmosphere truth field, since MERRA-2 stops at `64 km` and the
   residual above it is zero by construction.

## Resolved Julia Execution Blocker

Julia is available on the host machine, but this Codex shell has had repeated
filesystem-boundary issues executing it.

Observed attempts:

- `/home/space-falcon-1/.juliaup/bin/julia`
  Result: permission denied from this shell.
- `/tmp/julia`
  Result: not visible in this shell's `/tmp`.
- `/var/lib/snapd/hostfs/tmp/julia`
  Result: visible to `ls/stat`, but not readable or executable from this shell.
- `/home/space-falcon-1/Documents/SpaceAGORA.jl/tmp-julia`
  Result at time of check: this path existed as a directory, not an executable file.

The last concrete unblock requested was to place a runnable Julia executable in
the workspace at a file path such as:

- `/home/space-falcon-1/Documents/SpaceAGORA.jl/tmp-julia-bin`

## User-Provided Julia Paths During Debugging

The user indicated these Julia-related paths during the conversation:

- `/home/space-falcon-1/.juliaup/bin/julia`
- `/tmp/julia`
- workspace-local location under the repository root

Julia 1.12.1 is now usable from this shell.

## WamIPE

WamIPE's public density grid covers 100-1000 km and its installed 0.1
implementation omits internal download and NetCDF-cache helpers referenced by
its density API. MERRA-2 now supplies the lower atmosphere to `64 km`, so WamIPE
is the natural candidate for the `64-125 km` band the MERRA-2 truth source
leaves as a zero residual — but only if the implementation gap is closed.
