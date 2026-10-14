# Gravity harmonics inputs

CSV files with fully normalized spherical-harmonic coefficients, columns
`degree,order,C,S`, read by the harmonics gravity model through the
`gravity_harmonics_file` scenario key. Degree 0 and 1 rows are zero (the
point-mass term is applied separately).

| File | Body | Source |
|---|---|---|
| `EarthGGM05C.csv` | Earth | GGM05C |
| `egm96.csv` | Earth | EGM96 |
| `EGM96_GMAT_L50.csv` | Earth | EGM96 as distributed with GMAT R2025a (`data/gravity/earth/EGM96.cof` in the release tarball), degree/order 50, with that file's GM and radius header. Selected for the reported GMAT comparisons (`test/gmat_scenario_matrix.jl`, target `:gmat`). Coefficients are identical to `egm96.csv`. |
| `EGM96_GMAT_L50_zerotide.csv` | Earth | `EGM96_GMAT_L50.csv` with only C(2,0) converted to the zero-tide convention by the IERS 2010 permanent-tide term (−4.200675e-9). Selected for the STK comparison (target `:stk`); the retained reference runs' tide convention is unverified. |
| `Mars50c.csv` | Mars | Mars50c |
| `GMM2B.csv` | Mars | GMM-2B |
| `MGNP180U.csv` | Venus | MGNP180U |
| `LP165P.csv` | Moon | LP165P |
| `LP165P_permtide.csv` | Moon | `LP165P.csv` with the Earth-raised permanent tide added to C(2,0) and C(2,2) (IERS 2010 form; k2 = 0.02405, GRAIL GL0660B; not fitted). Selected for the STK comparison (target `:stk`); this does not establish the reference runs' tide convention. |
| `titan5.csv` | Titan | NASA GSFC PGDA product 91, `sha.titan_unnormalized` (Goossens et al. 2024, Nature Astronomy, doi:10.1038/s41550-024-02253-4), degree and order 5, zero-tide convention. Converted from unnormalized to fully normalized coefficients on 2026-09-04; header values GM = 8.978127e+12 m^3/s^2, reference radius = 2575000 m. |

`titan5.csv` is the file the aerobraking-perturbation study
(`benchmarks/studies/aerobraking_perturbation_mc`) names for Titan; the file
the study was originally run with was never tracked, so results with this
field may differ from the published ones.

The coefficient identities and the reference trajectories' generation settings
are separate evidence. The available older GMAT generator selects JGM2 for
Earth and has no Luna central-body case, while this comparison selects EGM96
and LP165P. The referenced full-arc rerun record has not been recovered. Until
scripts, input files and tool settings are pinned to the actual reference
bytes, GMAT/STK results are conditional comparisons. Reduced residuals alone
cannot establish a gravity or tide convention, nor isolate the entire error to
body-fixed rotations. Diagnostic frame substitutions are sensitivity results.
