# Implementation plan: 4,000-satellite interlink scheduling

Date: 2026-10-03

## Professor's request and acceptance criteria

The request is a scaling experiment and a targeted optimization, not a change to
the physical objective.

- **REQ-001:** Generate a seeded randomized Walker-Delta constellation with
  exactly 4,000 satellites, obtain a nonempty valid schedule, and save timings,
  allocation measurements and schedule output.
- **REQ-002:** Reuse array storage for availability and signed availability
  transitions. Preserve range, energy, temperature, activity and custom
  eligibility semantics.
- **REQ-003:** Avoid dictionary-heavy graph rebuilding and exponential matching
  on tree graphs. Preserve exact penalized matching, terminal exclusivity,
  deterministic selection, accepted-step scheduling and continuous forces.

## Technical context and constraints

Julia 1.12; existing SpaceAGORA, OrdinaryDiffEq, StaticArrays, CSV and DataFrames
dependencies. No new dependencies or network services are needed.

Source evidence:

- [Interlink model](I_src/1.4_satellite_interlink/inter_link_models.jl):
  registration and mutable connection state currently live in a dictionary.
- [Scheduler](I_src/1.4_satellite_interlink/scheduling_policies.jl):
  exact matching currently memoizes subsets of the whole terminal graph.
- [Force evaluation](I_src/1.3_dynamics/coupled/force_torque_models/laser_force_effectors.jl):
  every spacecraft currently scans all candidate edges on every RHS evaluation.
- [Existing example](II_examples/gve_sma_interlinks.jl) and
  [regressions](test/runtests.jl) define the physics and callback behavior.

No separate feature specification, constitution, topology or migration design
artifacts were supplied. This document records the user requirements and
source-based plan for this bounded Julia optimization, rather than inventing a
migration architecture. No framework migration guidelines apply.

Preserve the dictionary lookup API and mutable connection state. Topology must
be built through `register_candidate!`; cached arrays reference the same
connection objects. Keep all changes inside this copy. Do not introduce a
silent greedy approximation.

## Implementation steps and tasks

### P1.1: Reusable edge storage and transitions (REQ-002, REQ-003)

Use an edge-indexed array rather than a dense 4,000 by 4,000 matrix. Only
registered edges need storage. Keep per-satellite incident edge lists and
preallocate availability and signed transition arrays. A transition is
`Int8(new_available) - Int8(old_available)`; it is not a force or time evolution
matrix. Evaluate geometry without allocating a position-difference vector.

- [x] T001 [Plan:1.1] Add edge arrays, adjacency and bulk availability refresh in
  `I_src/1.4_satellite_interlink/inter_link_models.jl`.
- [x] T002 [Plan:1.1] Use incident edges in
  `I_src/1.3_dynamics/coupled/force_torque_models/laser_force_effectors.jl`.

### P2.1: Exact component matching (REQ-003)

Cache integer terminal endpoints and adjacency after registration. Reuse
iterative DFS, parent, order and dynamic-programming arrays. Split the
positive-net-score graph into components. Solve trees exactly in linear time
using bottom-up matching DP and iterative reconstruction. Retain exact subset
matching for cyclic components, using integer indices rather than keyed edge
lookups. Its exponential worst case remains explicitly documented.

- [x] T003 [Plan:2.1] Implement component/tree matching and array scoring in
  `I_src/1.4_satellite_interlink/scheduling_policies.jl`.
- [x] T004 [Plan:1.1,2.1] Add focused regression and exhaustive small-graph
  comparisons in `test/scalability.jl`; include them in `test/runtests.jl`.

### P3.1: Walker experiment and verification (REQ-001, REQ-002, REQ-003)

Use Walker-Delta T/P/F = 4000/80/F with seeded random phasing parameter,
common RAAN offset and common anomaly offset; preserve regular Walker spacing.
Use a circular 550 km orbit at 53 degrees and one terminal per satellite.
Register all 3,999 target-partner links: other edges always score zero under
the existing target-only `gve_sma` objective. Use a documented 2,000 km range
for this experiment (not a change to the model's 200 km default).

First measure scheduling alone over analytic two-body snapshots, with no claim
that this substitutes for force-coupled integration. Then run a short full-state
Tsit5 integration of all 4,000 spacecraft with actual callbacks and forces.
Keep output cadence and duration bounded to avoid large trajectory files.

- [x] T005 [Plan:3.1] Add seeded generator, snapshot benchmark, schedule export
  and full-engine smoke run in `II_examples/walker_interlinks.jl`.
- [x] T006 [Plan:3.1] Run existing physics/callback regressions and the 4,000
  satellite experiment; document measured results and limitations in
  `scalability-validation.md`.

## Validation strategy

- Compare exact total scores with exhaustive enumeration on small trees and
  cyclic/disconnected graphs; assert exclusivity and deterministic selection.
- Check range boundaries, availability transitions, activity and battery
  changes, cache reuse, and registration after a previous schedule.
- Run a 4,000-node star, chain and disconnected matching workload to verify
  there is no deep tree recursion.
- Compare incident-edge forces to the original all-edge definition.
- Run the existing suite including the one-hour physical step-size comparison.
- Warm up before measuring time and allocation; report observations, not
  hardware-independent latency guarantees.
- Assert exactly 4,000 spacecraft, nonempty schedules, finite state and
  diagnostics, terminal exclusivity, final integration time and one history
  entry per accepted step plus initialization.

## Requirement mapping

| REQ ID | Description | Plan items | Implementation evidence |
|---|---|---|---|
| REQ-001 | Reproducible 4,000-satellite schedule | P3.1 | Walker example, exported schedule and validation report |
| REQ-002 | Preallocated availability transitions | P1.1, P3.1 | Edge workspace and transition/reuse tests |
| REQ-003 | Faster exact graph updates and forces | P1.1, P2.1, P3.1 | DFS/tree solver, incident-edge forces, regression and scaling results |
