# Gravity harmonics inputs

CSV files of fully normalized spherical-harmonic coefficients, with columns
`degree,order,C,S` (`egm96.csv` adds `C_error,S_error` and `MGNP180U.csv` adds
`C_uncertainty,S_uncertainty`). `GravitationalHarmonicsModel(degree, order,
file, planet)` reads them; the verification study's scenarios name a file with
the `gravity_harmonics_file` key. The degree 0 and 1 rows are zero, because the
point-mass term is applied separately.

Most files start with a `# gravity_harmonics_source_metadata_v1` comment block
that records the source, citation, reference radius and GM, and the model uses
that radius and GM. `egm96.csv` and `titan5.csv` have no such block: unless the
caller passes `reference_radius_m` and `gm_m3s2`, the model falls back to the
planet's radius and GM, which may not be the values the coefficients were
solved with.

| File | Body | Source |
|---|---|---|
| `EarthGGM05C.csv` | Earth | GGM05C |
| `egm96.csv` | Earth | EGM96 |
| `EGM96_GMAT_L50.csv` | Earth | EGM96 as distributed with GMAT R2025a (`data/gravity/earth/EGM96.cof` in the release tarball), degree and order 50, with that file's GM and radius. Its coefficients equal those of `egm96.csv`. Used for the GMAT comparisons in `test/gmat_scenario_matrix.jl` (target `:gmat`). |
| `EGM96_GMAT_L50_zerotide.csv` | Earth | `EGM96_GMAT_L50.csv` with only C(2,0) moved to the zero-tide convention by the IERS 2010 permanent-tide term (−4.200675e-9). Used for the STK comparison (target `:stk`); the tide convention of the STK reference runs is not known. |
| `Mars50c.csv` | Mars | Mars50c |
| `GMM2B.csv` | Mars | GMM-2B |
| `MGNP180U.csv` | Venus | MGNP180U |
| `LP165P.csv` | Moon | LP165P |
| `LP165P_permtide.csv` | Moon | `LP165P.csv` with the Earth-raised permanent tide added to C(2,0) and C(2,2) (IERS 2010 form; k2 = 0.02405 from GRAIL GL0660B; not fitted). Used for the STK comparison (target `:stk`); that does not establish the tide convention of the reference runs. |
| `titan5.csv` | Titan | NASA GSFC PGDA product 91, `sha.titan_unnormalized` (Goossens et al. 2024, Nature Astronomy, doi:10.1038/s41550-024-02253-4), degree and order 5, zero-tide convention, converted from unnormalized to fully normalized coefficients. The file carries no GM or radius; the product's values are GM = 8.978127e+12 m^3/s^2 and reference radius 2575000 m, which `scripts/dev/viewer_demos/cassini_titan_flyby.jl` passes explicitly. |

`titan5.csv` is the Titan file the aerobraking-perturbation study (now in the
SpaceAGORAAtmosphereStudies research repository) names. The file that study
was first run with was never tracked, so results with this field may differ
from the published ones.

The GMAT and STK comparisons depend on more than these coefficients. The older
GMAT generator scripts that are available select JGM2 for Earth and have no
Moon-centred case, while these comparisons select EGM96 and LP165P, and the
record of the full-arc reference rerun has not been found. Until the scripts,
input files and tool settings are pinned to the exact reference files, treat
the GMAT and STK results as conditional comparisons. Smaller residuals alone
cannot establish a gravity or tide convention, nor attribute the whole error to
body-fixed rotations; runs with substituted frames are sensitivity results.
