# 5-Case Comparison

Requested case 3 duplicates case 2 and is run only once. Each case has one target and helpers at 1000 km with zero inclination. Both models use the same scheduler and 50 initial target orbital periods.

Trajectory recording interval: 1 second plus the final endpoint. Animation: disabled. Simulation settings otherwise inherit TEST_CASE from run_comparison.jl.

Wall times include model setup, simulation, recording, analysis, and any requested rendering; worker startup/imports are excluded, first-call compilation is included.

## case1_N20_h1050_i0_gve_sma

Target: 1050.0 km, 0.0 deg; helpers: 20; scheduler: gve_sma.

[Comparison report](comparison_v2_case1_N20_h1050_i0_gve_sma.md)

| Model | Simulated duration (s) | Model run wall time (s) |
| --- | --- | --- |
| ORACLE prototype | 318567.042 | 106.418467 |
| SpaceAGORA | 318567.042 | 121.094978 |

## case2_N20_h1000_i5_gve_sma

Target: 1000.0 km, 5.0 deg; helpers: 20; scheduler: gve_sma.

[Comparison report](comparison_v2_case2_N20_h1000_i5_gve_sma.md)

| Model | Simulated duration (s) | Model run wall time (s) |
| --- | --- | --- |
| ORACLE prototype | 315355.97 | 110.825029 |
| SpaceAGORA | 315355.97 | 117.452003 |

