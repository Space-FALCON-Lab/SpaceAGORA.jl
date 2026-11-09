# Kuang: later-task results

All assertions passed. Reviewer: Nihal. Generated: 2026-11-08T11:59:35.128. Runtime: 44.07 s.

## Shortened forty-spacecraft momentum run

Adapted from the existing forty-spacecraft example: 30 helpers at 1000 km, 10 targets at 1100 km, circular equatorial orbits, original ring spacing and 1 degree target phase offset. All masses remain 227 kg and each spacecraft has one virtual terminal. Allowed pairs are all 300 helper-to-target pairs; other pairs are not registered. This run does **not** replace the earlier unequal-mass check.

- Duration: 1608.944572739 s (0.25 target orbits, rather than 20).
- Recorded scheduling samples: 168; selected-set changes: 1; maximum active links: 10.
- Target SMA change relative to a separately propagated no-link control: 12.8810229 to 12.881023 m.
- Maximum no-link SMA drift across all 40 spacecraft: 6.66826963e-07 m (limit 0.02 m).

`LaserMomentumLink`: 10,000 W, magnification 100, range 200 km; ideal equal/opposite force `100 * 10000 / 299792458` N. `SemimajorAxisRatePolicy(targets=31:40)` is checked against the independent energy-rate formula `2*a^2*dot(v,F)/(mu*m)` at every recorded state. All feasible terminal components were verified to be disjoint edges, so selecting all positive-score edges is the exact optimum even though the full allowed graph is cyclic. No Python backend was loaded.

Attitude, animation, plotting and atmosphere are disabled. The ideal momentum model uses range-only feasibility: no mounted-beam capture, pointing or occultation is claimed here. Gaussian capture and pointing loss exist in the separate beam model; their ORACLE adequacy belongs with Nihal's physics work.

Evidence: [settings](results/constellation/settings.toml), [pre-run candidates](results/constellation/candidates_before_run.csv), [every candidate at every accepted step](results/constellation/candidate_scores.csv), [matching checks](results/constellation/accepted_steps.csv), [recorded contact/active windows](results/constellation/link_schedule.csv), [SMA comparison](results/constellation/spacecraft_comparison.csv), [linked states](results/constellation/states.csv), [no-link states](results/constellation/control_states.csv).

## Separate data-link experiments

Circular co-orbital spacecraft at 700 km, separated by 10 km of arc. All endpoints have finite data stores; the final spacecraft is a finite-store **destination**, not a ground station. Directed `DataLink`s use `FixedRateBudget`, virtual terminals and zero transmit power to isolate scheduling from RF/power physics. Data-link planet occultation is enabled by the built-in model and independently checked. All runs last 20 s with backlog horizon 5 s.

| Case | dtmax (s) | Samples | Set changes | Delivered to destination (bit) | Production (bit) | Permitted production overflow (bit) | Max node balance residual (bit) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [producer_control](results/backlog/producer_control/accepted_steps.csv) | 0.250 | 85 | 0 | 4000000.000000000 | 4000000.000000000 | 0.000000000 | 0 |
| [production_competition](results/backlog/production_competition/accepted_steps.csv) | 0.250 | 87 | 9 | 11082065.734343071 | 4000000.000000000 | 0.000000000 | 9.31323e-10 |
| [busy_relay](results/backlog/busy_relay/accepted_steps.csv) | 0.250 | 87 | 1 | 20000060.650488533 | 20000000.000000000 | 19999878.699022926 | 1.49012e-08 |
| [direct_route_control](results/backlog/direct_route_control/accepted_steps.csv) | 0.250 | 86 | 0 | 29999999.999999996 | 19999999.999999996 | 0.000000000 | 1.49012e-08 |
| [production_competition_refined](results/backlog/production_competition_refined/accepted_steps.csv) | 0.125 | 167 | 19 | 11082065.734343072 | 3999999.999999999 | 0.000000000 | 9.31323e-10 |
| [busy_relay_refined](results/backlog/busy_relay_refined/accepted_steps.csv) | 0.125 | 166 | 1 | 20000060.650488533 | 20000000.000000000 | 19999878.699022926 | 1.49012e-08 |

### What the fixtures establish

1. **Empty producer control:** a 200,000 bit/s source remains selected and delivers exactly 4,000,000 bit over 20 s, to tolerance.
2. **Competition:** that same empty producer initially has a positive score, but loses the destination's single terminal to another sender with 8,000,000 bit. Later selections and residual backlog are recorded: ongoing production is not a selection guarantee, and backlog score is not a fairness guarantee.
3. **Busy relay:** source 1 starts with 80,000,000 bit; relay 2 starts full at 20,000,000 bit and produces 1,000,000 bit/s; destination 3 starts empty. Rates are 1,000,000 bit/s for source-to-relay and relay-to-destination, and 500,000 bit/s direct. Terminal counts are (1,2,2). The terminal graph is a forest: relay input and output use distinct terminals, as do the destination's direct and relay inputs.

Initially the full relay has no **previous** outflow credit, so the scheduler selects direct + relay-out. At the next accepted step, its previously active feasible relay-out supplies 1,000,000 bit/s credit; relay-in + relay-out scores 2,000,000 bit/s, beating direct + relay-out at 1,500,000 bit/s. Actual delivery then falls to 1,000,000 bit/s, while the relay's own production overflows at 1,000,000 bit/s. Incoming traffic displaces local production; a full receiver not overfilling is not a no-loss proof.

The same-backlog-policy direct-route control explicitly disables relay-in, delivering 9999939.349511463 bit more and losing 19999878.699022926 bit less than the busy-relay run. This is a feasible routing comparator, **not** a claim of a globally optimal routing policy. Store counters do not track packet origin, so destination totals must not be described as source-1-only delivery.

Scores are checked independently as `min(rate, held_tx/5 + production_tx, room_rx/5 + previous_feasible_outflow_rx)`. All matching subsets are enumerated at every scheduling sample. A read-only scorer replay also checks rejected candidates, using the previous selection and store flags reconstructed at that state, never the newly selected links. The separate [outflow probes](results/backlog/outflow_probes.csv) demonstrate zero credit with no previous link, positive credit with a previous feasible link, and zero credit when that previous link is out of range.

Every case saves its spacecraft/stores, pre-run candidates and matching subsets, full candidate-score audit, selected-set history, contact/active windows, per-node accounting and network accounting. [Analytic comparisons](results/backlog/hand_comparisons.csv) integrate the simple constant-rate modes over the **recorded accepted-step** intervals, not a fictional continuously located scheduling event. [Step-refinement results](results/backlog/step_refinement.csv) compare dtmax 0.25 s with 0.125 s for competition and relay selection. Numerical solver error and scheduler time-resolution sensitivity are different quantities.

## Accounting, frames and tolerances

- Positions/velocities are planet-centered inertial Cartesian SI values; initial angles are degrees. No attitude or body torques are propagated.
- Solver: Tsit5; reltol 1e-10, abstol 1e-8; dtmax 10 s for the constellation, 0.25 s for data cases, 0.125 s for refined cases. Every accepted state is retained.
- Momentum scores: rtol 1e-12, atol 1e-14 m/s. Data scores: rtol 1e-12, atol 1e-5 bit/s. Matching ties may select any equal-optimal set.
- Data counter/hand-comparison absolute tolerance: 1 bit. Store bounds and monotonic counters are checked at every saved accepted state. Store scoring uses the engine's documented full/empty threshold of 1e-9 times capacity.
- Per node: `stored = initial + produced + received - sent - lost`. Internal sent and received totals must agree. Network accounting includes the destination store, so delivery is **not subtracted again**: `sum(initial) + sum(produced) = sum(stored) + sum(lost)`.
- Loss is permitted only at relay 2 in the busy-relay cases, with its quantity independently predicted. Other fixtures require zero loss to 1 bit. These bounded checks make no global no-loss claim.
- `record_available=true`; the saved schedule is the engine's history converted by `link_schedule`, not a schedule inferred from raw force/rate queries. Windows and switches have accepted-step resolution.

## Coordination and separate work

Give Nihal the per-node/network accounting and hand-comparison tables for review of production, transfers and permitted overflow; Kuang owns the matching and delivery comparison. Ground-station elevation/planet geometry, viewer/animation validation, a genuinely cyclic feasible-graph backend test, mounted-beam physics and resumed laser history are **not tested** here. Resumed laser history remains Nihal's later task. Fixed rates, no acquisition time and no hardware power budget are deliberate scheduler fixtures, not a validated ORACLE communication model.

## Reproducibility

Julia 1.12.1; SpaceAGORA 0.7.1; checkout `70e6cfde4efeab9a8129b523dde799f7c723bd1e`. [Provenance and hashes](results/provenance.toml), [run status](results/status.toml), [data summary](results/backlog/summary.csv), [run instructions](README.md). No existing repository files were edited and no dependencies were installed. All generated deliverables are inside this folder.
