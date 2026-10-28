# Bounded capability verdicts

See REPORT.md for exact scientific sources, the dispatch-equivalent follow-up
KS repair and the final installation pin. Passing CI alone is not these verdicts.

| Capability | Verdict | Limit/evidence |
|---|---|---|
| One-spacecraft two-body24h | Fine PASS; coarse/medium FAIL |Analytical Kepler plus independent refined Cartesian reference;1,441 epochs|
| Static-axis J2-only24h | Fine PASS; coarse/medium FAIL |Specified SI spherical datum,+z symmetry axis|
| Co-rotating exponential drag24h | Fine PASS; coarse/medium FAIL |Constant coefficients, no wind/history/frozen sample semantics|
| J2 plus drag24h | Fine PASS; coarse/medium FAIL |Same frozen physical inputs and mass/heat gates|
| Full configured acceleration and auxiliary coupling | PASS tested states |68 parity checks; mass/heat in same nonlinear stage, one monopole|
| Physical-time control/actions/ideal burn/impulse | PASS specified6000s case |13s control,17s output; independent Cartesian recomputation|
| Continuous observer root and retained dense history | PASS tested patterns |Root refinement; command/impulse/knot/derivative/rejection checks; no arbitrary-event guarantee|
| Segmented and fresh-process resume | PASS specified cases |Full persisted prefix; no missing/duplicate nominal epochs; state/history restored|
| Coincident221s action/control/output/checkpoint | PASS442s case |Impulse then control then output, no duplicate committed action|
| Normal installation with/without KS | PASS checkout/installer route |Fresh maintained-repository source, actual assets;316/9 runtime checks|
| Direct Pkg.add SpaceAGORA Git URL | FAIL observed installation route |LibGit2 submodule checkout error; use verified checkout/installer|
| Unsupported configuration refusal | PASS selected guards |Includes new data owners/effective central scope; custom purity trait is a contract|
| Fixed-physical-time sensitivities/tangent map | Not implemented/claimed |Standalone derivative results do not qualify this new reset map|
| Integrated sails/high-order gravity/multiple spacecraft/interlinks | Deferred |No integration verdict inherited from standalone examples|
| Attitude/articulation/mounted systems/history-mutating RHS | Outside initial support |Relevant configurations refused|
| TrajectoryRecorder cadence/specialized split-multirate output | Outside table-output qualification |Separate work required|
| MPC predictor/optimizer/feedback changes | Unchanged |No replacement or integration|
| Two-body matched-accuracy timing | Bounded observed result |About106s KS vs4.81ms DP8;3zero-compile warm calls; machine/noise limits|
| J2-plus-drag matched-accuracy timing | Inconclusive wall performance |Accuracy passes; ten long monitor gaps; ratio withheld, all calls retained|

Fine: relative1e-12, absolute1e-13, maximum physical step5s for these inputs.
Targets: maximum sampled1m/0.001m/s plus declared mass/heat/clock budgets. All
eight coarse/medium misses remain in evidence/baseline-results.toml. Broader
support needs separate declared verification.
