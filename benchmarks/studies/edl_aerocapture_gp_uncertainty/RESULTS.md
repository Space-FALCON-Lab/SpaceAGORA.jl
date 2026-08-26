# Results

Generated on 2026-08-23T21:40:34 UTC.

Truth source: `merra2_native` (position-indexed: true).
Prior source: `gram`.

> Prior and truth share a lineage: GRAM's nominal below about `65 km` is
> the MERRA-2 climatology, so the residual here is a specific day minus
> its own monthly mean. Run `--prior-source nrlmsise00` for an
> independent prior.
Aerocapture sensing band: `75.0` to `35.0 km`.
GRAM correlation scales at `35.0 km`: horizontal `2279.4 km`, vertical `18.33 km`, temporal `26.64 hr`.
Kernels: `squared_exponential`, `matern32`, `matern52`, `gram_exponential`.
Mean bases: `none`, `constant`, `linear_alt`.
Amplitude mode: `prior_scaled`, lambda bounds `auto`.

Means over 36 cases. Nominal coverage is `0.683` and `0.954`.

## GRAM Baseline

| weighted log RMSE | weighted NLPD | coverage `1-sigma` | coverage `2-sigma` |
|---|---|---|---|
| 0.033976 | -5.3073 | 0.35 | 0.405 |

The baseline is well under-dispersed: GRAM's ensemble is roughly half as
wide as the truth field requires. Every model below inherits that.

## Models, Ranked by Mean Weighted Log-Density RMSE

| parameterization | kernel | mean basis | vs GRAM | weighted NLPD | cov `1-sigma` | cov `2-sigma` |
|---|---|---|---|---|---|---|
| `density_scale_factor` | `gram_exponential` | `constant` | 45.64% | -8.24 | 0.536 | 0.739 |
| `log_density_scale_factor` | `gram_exponential` | `constant` | 44.58% | -8.73 | 0.568 | 0.783 |
| `log_density` | `gram_exponential` | `constant` | 44.58% | -8.73 | 0.574 | 0.783 |
| `density_scale_factor` | `matern32` | `constant` | 38.73% | -8.12 | 0.505 | 0.717 |
| `density_scale_factor` | `matern52` | `constant` | 38.43% | -8.07 | 0.496 | 0.7 |
| `log_density` | `matern32` | `constant` | 38.35% | -8.63 | 0.535 | 0.751 |
| `log_density_scale_factor` | `matern32` | `constant` | 38.35% | -8.63 | 0.53 | 0.751 |
| `density_scale_factor` | `gram_exponential` | `none` | 35.57% | -8.25 | 0.477 | 0.696 |
| `density_scale_factor` | `matern32` | `none` | 33.8% | -7.91 | 0.515 | 0.719 |
| `log_density_scale_factor` | `matern52` | `constant` | 33.2% | -8.59 | 0.522 | 0.74 |
| `log_density` | `matern52` | `constant` | 33.2% | -8.59 | 0.527 | 0.74 |
| `log_density_scale_factor` | `gram_exponential` | `none` | 31.91% | -8.58 | 0.521 | 0.758 |

Accuracy and calibration do not rank together: read both columns before
picking a configuration.

## By Mean Basis

| mean basis | vs GRAM | weighted NLPD | cov `1-sigma` | cov `2-sigma` |
|---|---|---|---|---|
| `none` | -0.69% | -7.89 | 0.527 | 0.742 |
| `constant` | 34.85% | -8.45 | 0.523 | 0.736 |
| `linear_alt` | 8.97% | -8.45 | 0.445 | 0.756 |

## By Aerocapture-to-EDL Gap

| gap | vs GRAM |
|---|---|
| 1.0 hr | 20.16% |
| 3.0 hr | 5.43% |
| 6.0 hr | 17.51% |

## Fitted Mean Function

`beta_constant` is the bulk log-density offset the aerocapture pass
inferred relative to the GRAM prior. Unlike the zero-mean GP term it
extrapolates into the unsensed band.

| parameterization | mean basis | beta_constant | beta_alt_slope | lambda | signal fraction |
|---|---|---|---|---|---|
| `log_density` | `constant` | 0.03011 | 0.0 | 0.8434 | 0.14 |
| `log_density` | `linear_alt` | 0.02487 | -0.01823 | 0.6116 (some at bound) | 0.095 |
| `density_scale_factor` | `constant` | 0.02608 | 0.0 | 0.877 | 0.137 |
| `density_scale_factor` | `linear_alt` | 0.02091 | -0.01853 | 0.6476 (some at bound) | 0.091 |
| `log_density_scale_factor` | `constant` | 0.03011 | 0.0 | 0.8434 | 0.14 |
| `log_density_scale_factor` | `linear_alt` | 0.02487 | -0.01823 | 0.6116 (some at bound) | 0.095 |

`signal fraction` is the share of training-residual variance above the
`5%` measurement-noise floor. Values near zero mean the pass carried
little information and `lambda` is barely identifiable.

## Output Files

- `output/edl_gp_e75/summary_metrics.csv`
- `output/edl_gp_e75/kernel_comparison.csv`
- `output/edl_gp_e75/mean_function_fits.csv`
- `output/edl_gp_e75/cases.csv`
- `output/edl_gp_e75/pointwise_predictions.csv` when `--save-pointwise=true`
