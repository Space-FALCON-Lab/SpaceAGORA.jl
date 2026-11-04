# Absolute R/N campaign validation

Checked UTC: 2026-10-22T00:47:02.133.

## Scope and source preservation

- Exactly **270 fresh cases**, 100 periods each; 1001 rows with endpoints.
- Targets 900/950/1000/1050/1100 km; inclinations 0/0.5/1 degree.
- Helpers initially at 1000 km; counts 1/100/200/300/400/500/600/700/800; both schedulers.
- Row 801 is the 80-orbit endpoint. No separate 80-orbit runs.
- Original runner files match the pre-change snapshot. Old N700 remains stopped, 4 verified cases preserved.
- [Source archive](source_snapshot_pre_absolute_RN/sources.tar.gz): 866 checksummed files; SHA256 `bcebb39026f30a770d22cfa8d31d1b67ec92a4660ad99f1efb56b6674f045ab8`.
- Shared library edits are limited to diagnostic RHS/initialization, default saving and checkpoint validation.

## Tests and independent validation

- **9085 assertions passed** in [the focused suite](fixtures_final.log), including existing interlink/scheduler tests and enabled checkpoint round trips.
- **13 pilots/baselines** passed fresh-file verification; pilot data remain in [PilotOutput](../PilotOutput).
- Net and absolute integrals are independently reconstructed from dense Cartesian states and the actually held link, not saved diagnostics.
- Absolute acceleration is integrated directly, with zero-crossing splitting and 16-versus-8-panel Simpson convergence checks.
- Maximum independent net error: 7.7e-14 m/s; absolute-integral error: 8.02e-09 m/s; summed quadrature difference: 5.26e-13 m/s.
- R/N branch tolerances: 1e-8 m/s; existing T allowance: 1e-7 m/s; net-sum tolerance: 1e-8 m/s. No values were clamped and no limits were relaxed.
- Independent quadrature covers the short instrumented pilots. Their measured errors are not claimed as a uniform 100-orbit cumulative-error bound; production still checks every saved row and the physical/event invariants.
- Inclined lower/higher-target cases exercise both R signs for each scheduler; N has both signs. See [sign coverage](sign_coverage.toml).
- Zero-power and coplanar inclination pilots have exactly zero laser diagnostics. Recorder/quadrature instrumentation leaves the same-layout dynamics and schedule bitwise unchanged.

## Frozen-library physical comparisons

| Comparison | Maximum position difference (m) | Maximum velocity difference (m/s) | Held-schedule disagreement (s) |
|---|---:|---:|---:|
| inc1050_inc1050_control | 0 | 0 | 0 |
| baseline_sma800_new_sma800 | 96.5949 | 0.0965517 | 196.752 |
| baseline_inc800_new_inc800 | 71.15 | 0.0708912 | 183.269 |
| baseline_zero_zero_sma | 5.3743e-08 | 8.12714e-11 | 0 |

Added states may change adaptive steps. Let a be maximum laser acceleration, D the measured time when held links differ, T the duration and n the initial mean motion. The acceleration difference is at most 2a during those intervals. Conservative near-circular comparison tolerances are 0.02 + 10*T*a*D metres and 5e-5 + 14*a*D + 1.01*n*(10*T*a*D) m/s. The latter includes rotation of the velocity vector as phase changes. The schedule-disagreement limit is 2%. These are comparison allowances, not solver tolerances. The RHS tests separately require identical physical derivatives at identical states.

## Serial 800-helper benchmarks

| Implementation/policy | Orbits | Simulation wall (s) | Post-warmup wall/orbit (s) | Peak RSS (GiB) | Temporary spool (MiB) |
|---|---:|---:|---:|---:|---:|
| baseline_sma800 | 10 | 424.20 | 7.98 | 2.188 | 231.38 |
| new_sma800 | 10 | 423.16 | 8.01 | 2.203 | 231.38 |
| baseline_inc800 | 10 | 426.51 | 8.30 | 2.201 | 231.38 |
| new_inc800 | 10 | 423.59 | 8.06 | 2.205 | 231.38 |

Updated campaign estimate: **40.6 hours**; planning allowance: **46-50 hours**.

The estimate starts with 240 observed old cases plus 30 N700 interpolations, adds the existing process/verification allowance to recover the 41-hour baseline, and applies each scheduler's measured old/new 100-orbit projection ratio. A projection counts startup once, then adds 90 periods at the post-warmup rate. It does not multiply a short pilot's entire wall time by ten or fifty. Hardware load and case-to-case differences remain uncertain.

## Capacity and launch readiness

- Available disk: 48.60 GiB. Expected persistent output: 4.95 GiB; conservative allowance: 10.47 GiB.
- Largest one-case spool: 2.31 GiB; total disk requirement including the 5 GiB reserve: 17.78 GiB.
- Largest measured pilot RSS: 2.21 GiB; available memory: 20.14 GiB; required margin: 6.00 GiB.
- No existing data was deleted to obtain capacity. In-run and pre-case disk guards remain enabled.
- Validated runtime source signature: `0834b044dd9183906c2e8058796f5cc70df5bff329000102ca8c8c53b1c7de48`.

Preflight passed. Actual launch and the first production verification are recorded separately in [launch confirmation](launch_confirmation.toml). Monitor/stop/resume commands are in the [runner README](../../../Simulation_absolute_RN/Code/README.md).
