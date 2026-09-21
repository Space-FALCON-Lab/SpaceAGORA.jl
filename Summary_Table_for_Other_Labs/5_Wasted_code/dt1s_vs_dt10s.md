# Saved dt1s versus dt10s Results

Existing outputs from both models are compared within each model, for cases 1 and 2. No dynamics were rerun and no trajectory was downsampled. Target altitude is labelled `t`; helper altitude is labelled `h`.

Checks cover complete output grids, all 21 spacecraft positions and velocities, shared timestamps including the final endpoint, active helpers, full recorded event tables, and encounter extrema. Endpoint tolerance is 1e-8 s; interior timestamps match exactly. Differences below are measurements, not acceptance tolerances.

Missing dt1s analysis CSVs are reconstructed from their original Feather bundles in temporary directories using the existing extractor. The dt10s analysis CSVs are read from the saved output folders. Original outputs are unchanged. Runtime and solver return codes cannot be recovered from these Feather bundles. Constants in the generated case summaries use the existing model configuration.

## case1_N20_t1050_i0_gve_sma

### prototype

[dt1s source Feather bundle](../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/)

[dt10s source Feather bundle](../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/)

| Cadence | Rows | Start (s) | End (s) | Grid and finite-state checks |
| --- | --- | --- | --- | --- |
| dt1s | 318569 | 0 | 318567.042 | PASS |
| dt10s | 31858 | 0 | 318567.042 | PASS |

| Quantity | Maximum absolute shared-time difference | Maximum absolute final difference |
| --- | --- | --- |
| Position components (m) | 0 | 0 |
| Velocity components (m/s) | 0 | 0 |
| Accumulated delta-v components (m/s) | 0.000140396055 | 5.67743064e-05 |

Shared timestamps: 31858; active-helper mismatches: 0.

| Event table | dt1s rows | dt10s rows | Exactly equal (all fields) |
| --- | --- | --- | --- |
| geometry_encounters.feather | 11 | 11 | true |
| laser_on.feather | 10 | 10 | true |

Encounter analysis tables exactly equal: true. Encounter-state rows: 53898 (dt1s), 5411 (dt10s).

| Metric | Extremum | dt1s | dt10s | dt10s minus dt1s |
| --- | --- | --- | --- | --- |
| relative_speed_mps | min | 16.2480891 | 16.249477 | 0.00138793103 |
| relative_speed_mps | max | 211.308821 | 211.308821 | 0 |
| range_rate_mps | min | -92.8854155 | -92.8854155 | 0 |
| range_rate_mps | max | 103.790955 | 103.790627 | -0.000328045247 |
| geometry_duration_s | min | 2591.86384 | 2591.86384 | 0 |
| geometry_duration_s | max | 5368.35988 | 5368.35988 | 0 |
| laser_on_s | min | 0 | 0 | 0 |
| laser_on_s | max | 3432.29205 | 3432.29205 | 0 |

### spaceagora

[dt1s source Feather bundle](../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/)

[dt10s source Feather bundle](../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/)

| Cadence | Rows | Start (s) | End (s) | Grid and finite-state checks |
| --- | --- | --- | --- | --- |
| dt1s | 318569 | 0 | 318567.042 | PASS |
| dt10s | 31858 | 0 | 318567.042 | PASS |

| Quantity | Maximum absolute shared-time difference | Maximum absolute final difference |
| --- | --- | --- |
| Position components (m) | 0 | 0 |
| Velocity components (m/s) | 0 | 0 |
| Accumulated delta-v components (m/s) | 0 | 0 |

Shared timestamps: 31858; active-helper mismatches: 0.

| Event table | dt1s rows | dt10s rows | Exactly equal (all fields) |
| --- | --- | --- | --- |
| geometry_encounters.feather | 11 | 11 | true |
| laser_on.feather | 11 | 11 | true |

Encounter analysis tables exactly equal: true. Encounter-state rows: 53894 (dt1s), 5411 (dt10s).

| Metric | Extremum | dt1s | dt10s | dt10s minus dt1s |
| --- | --- | --- | --- | --- |
| relative_speed_mps | min | 16.3945998 | 16.3965205 | 0.00192070927 |
| relative_speed_mps | max | 211.289426 | 211.289426 | 0 |
| range_rate_mps | min | -92.7039875 | -92.7039875 | 0 |
| range_rate_mps | max | 103.738403 | 103.738293 | -0.000110192474 |
| geometry_duration_s | min | 2591.86361 | 2591.86361 | 0 |
| geometry_duration_s | max | 5367.85861 | 5367.85861 | 0 |
| laser_on_s | min | 0.000379693605 | 0.000379693605 | 0 |
| laser_on_s | max | 3422.33646 | 3422.33646 | 0 |

[dt10s model comparison](comparison_v2_case1_N20_t1050_i0_gve_sma_dt10s.md)

## case2_N20_t1000_i5_gve_sma

### prototype

[dt1s source Feather bundle](../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1000km_ih0_it5_e0_nu0_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/)

[dt10s source Feather bundle](../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1000km_ih0_it5_e0_nu0_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/)

| Cadence | Rows | Start (s) | End (s) | Grid and finite-state checks |
| --- | --- | --- | --- | --- |
| dt1s | 315357 | 0 | 315355.97 | PASS |
| dt10s | 31537 | 0 | 315355.97 | PASS |

| Quantity | Maximum absolute shared-time difference | Maximum absolute final difference |
| --- | --- | --- |
| Position components (m) | 0 | 0 |
| Velocity components (m/s) | 0 | 0 |
| Accumulated delta-v components (m/s) | 0.000139401882 | 0.000102874063 |

Shared timestamps: 31537; active-helper mismatches: 0.

| Event table | dt1s rows | dt10s rows | Exactly equal (all fields) |
| --- | --- | --- | --- |
| geometry_encounters.feather | 101 | 101 | true |
| laser_on.feather | 20 | 20 | true |

Encounter analysis tables exactly equal: true. Encounter-state rows: 63402 (dt1s), 6523 (dt10s).

| Metric | Extremum | dt1s | dt10s | dt10s minus dt1s |
| --- | --- | --- | --- | --- |
| relative_speed_mps | min | 609.365442 | 609.365442 | 0 |
| relative_speed_mps | max | 643.92912 | 643.925453 | -0.00366709846 |
| range_rate_mps | min | -642.348908 | -642.326988 | 0.0219202347 |
| range_rate_mps | max | 642.349261 | 642.345065 | -0.00419596962 |
| geometry_duration_s | min | 317.175254 | 317.175254 | 0 |
| geometry_duration_s | max | 634.329102 | 634.329102 | 0 |
| laser_on_s | min | 0 | 0 | 0 |
| laser_on_s | max | 317.175254 | 317.175254 | 0 |

### spaceagora

[dt1s source Feather bundle](../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1000km_ih0_it5_e0_nu0_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/)

[dt10s source Feather bundle](../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1000km_ih0_it5_e0_nu0_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/)

| Cadence | Rows | Start (s) | End (s) | Grid and finite-state checks |
| --- | --- | --- | --- | --- |
| dt1s | 315357 | 0 | 315355.97 | PASS |
| dt10s | 31537 | 0 | 315355.97 | PASS |

| Quantity | Maximum absolute shared-time difference | Maximum absolute final difference |
| --- | --- | --- |
| Position components (m) | 0 | 0 |
| Velocity components (m/s) | 0 | 0 |
| Accumulated delta-v components (m/s) | 0 | 0 |

Shared timestamps: 31537; active-helper mismatches: 0.

| Event table | dt1s rows | dt10s rows | Exactly equal (all fields) |
| --- | --- | --- | --- |
| geometry_encounters.feather | 101 | 101 | true |
| laser_on.feather | 20 | 20 | true |

Encounter analysis tables exactly equal: true. Encounter-state rows: 63401 (dt1s), 6523 (dt10s).

| Metric | Extremum | dt1s | dt10s | dt10s minus dt1s |
| --- | --- | --- | --- | --- |
| relative_speed_mps | min | 609.36543 | 609.36543 | 0 |
| relative_speed_mps | max | 643.933606 | 643.929998 | -0.0036081526 |
| range_rate_mps | min | -642.348924 | -642.327001 | 0.021923509 |
| range_rate_mps | max | 642.349148 | 642.344893 | -0.00425552162 |
| geometry_duration_s | min | 317.175305 | 317.175305 | 0 |
| geometry_duration_s | max | 634.329225 | 634.329225 | 0 |
| laser_on_s | min | 0 | 0 | 0 |
| laser_on_s | max | 317.174925 | 317.174925 | 0 |

[dt10s model comparison](comparison_v2_case2_N20_t1000_i5_gve_sma_dt10s.md)

## Interpretation

Zero position/velocity differences mean identical saved dynamics at the shared times, not necessarily identical sampled extrema. A one-second grid can capture motion extrema missed by a ten-second grid. Prototype accumulated delta-v uses trapezoidal integration on saved samples, so its diagnostic can differ even when shared trajectory states and event intervals agree. This comparison does not claim the two different models are physically equivalent.

