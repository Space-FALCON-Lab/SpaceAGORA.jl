# Gravity harmonics inputs

CSV files with fully normalized spherical-harmonic coefficients, columns
`degree,order,C,S`, read by the harmonics gravity model through the
`gravity_harmonics_file` scenario key. Degree 0 and 1 rows are zero (the
point-mass term is applied separately).

| File | Body | Source |
|---|---|---|
| `EarthGGM05C.csv` | Earth | GGM05C |
| `egm96.csv` | Earth | EGM96 |
| `Mars50c.csv` | Mars | Mars50c |
| `GMM2B.csv` | Mars | GMM-2B |
| `MGNP180U.csv` | Venus | MGNP180U |
| `LP165P.csv` | Moon | LP165P |
| `titan5.csv` | Titan | NASA GSFC PGDA product 91, `sha.titan_unnormalized` (Goossens et al. 2024, Nature Astronomy, doi:10.1038/s41550-024-02253-4), degree and order 5, zero-tide convention. Converted from unnormalized to fully normalized coefficients on 2026-09-04; header values GM = 8.978127e+12 m^3/s^2, reference radius = 2575000 m. |

`titan5.csv` is the file the aerobraking-perturbation study
(`benchmarks/studies/aerobraking_perturbation_mc`) names for Titan; the file
the study was originally run with was never tracked, so results with this
field may differ from the published ones.
