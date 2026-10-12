# Walker constellation scaling: implementation and measured results

Date: 2026-10-03. Julia 1.12.1, Linux. Changes are confined to this copy.
See [plan.md](plan.md) for the professor's request, implementation tasks and
requirement mapping.

## What changed

1. [The model](I_src/1.4_satellite_interlink/inter_link_models.jl) now keeps
   edge-indexed keys and connection references, per-spacecraft incident edges,
   a reusable Boolean availability array and an `Int8` transition array.
   `update_availability!(model, spacecraft, state; is_active=...)` updates all
   candidates in place. The original single-key overload is retained.
2. [The scheduler](I_src/1.4_satellite_interlink/scheduling_policies.jl) caches
   a 2-by-E integer terminal-endpoint matrix, adjacency lists and traversal/DP
   arrays. It computes connected components of the available positive-net-score
   graph with iterative DFS. Trees use exact bottom-up matching and iterative
   reconstruction. Cyclic components retain exact subset matching, independently.
3. [Force evaluation](I_src/1.3_dynamics/coupled/force_torque_models/laser_force_effectors.jl)
   visits only the current spacecraft's incident edges, not the entire graph.
   This reduces an all-spacecraft pass from O(N*E) edge inspections to O(N+E).
   Stage geometry and the force law remain unchanged.
4. [The Walker example](II_examples/walker_interlinks.jl) provides a seeded
   generator, analytic snapshot benchmark, CSV exports and a full-engine smoke
   run with continuous interlink forces and accepted-step callbacks.

Use `register_candidate!` to change topology; do not insert, delete or replace
entries directly in `linkgraph` or edit workspace arrays. Dictionary lookups and
changes to existing connection **state** still reference the same objects used
by the arrays. Registration invalidates the matching workspace; the next
selection rebuilds it once.

For edge `e`, `available[e]` is the last refreshed value, and
`availability_transition[e]` is `new - old`: +1 entered availability, -1 left,
0 unchanged. These are overwritten at each refresh, not accumulated. The
single-key overload updates only that edge. Direct state edits become reflected
in these arrays on refresh. Availability does not depend on previous occupancy.

An edge array is preferable to a dense N-by-N availability matrix here: it
stores only configured candidate connections. Geometry/eligibility still require
O(E) evaluation; matrix addition alone cannot determine which links are in range.

## Reproduce

From the outer workspace root:

```sh
SPACEAGORA_RHS_CALIBRATE=off julia --project='8_SpaceAGORA.jl-main_v5 copy' --startup-file=no '8_SpaceAGORA.jl-main_v5 copy/test/runtests.jl'
SPACEAGORA_RHS_CALIBRATE=off julia --project='8_SpaceAGORA.jl-main_v5 copy' --startup-file=no '8_SpaceAGORA.jl-main_v5 copy/II_examples/walker_interlinks.jl'
```

The test command includes the existing one-hour small-system physics comparison.
The Walker command performs **both** a one-hour analytic snapshot experiment and
a separate 20-second full force-coupled integration. Snapshot propagation is
not a substitute for force-coupled propagation.

## Configuration

| Setting | Value |
|---|---|
| Walker-Delta T/P/F | 4000/80/25 |
| Satellites per plane | 50 |
| Random seed | 20261003 |
| Circular altitude / inclination | 550 km / 53 degrees |
| Random RAAN offset | 292.43964213039555 degrees |
| Random anomaly offset | 129.40018559361184 degrees |
| Target / terminals per satellite | Satellite vector index 1 / 1 |
| Registered candidates | All 3,999 target-partner connections |
| Example range | 2,000 km |
| Power / B | 10,000 W / 100 |
| Objective | Existing target-only `gve_sma`, zero active-link penalty |

Randomization changes the Walker phasing integer and common RAAN/anomaly offsets;
it does not destroy the regular Walker spacing. Model defaults, including the
original 200 km range, are unchanged. The longer range belongs only to this
experiment. With a single target terminal, **one selected link is the expected
optimal schedule**, not a failure to schedule all satellites. Non-target links
have zero score under this policy, so omitting them does not change the optimum.

## Results

### Snapshot scheduling: 4,000 satellites, 61 samples over one hour

| Measurement | Observed |
|---|---|
| Available links per snapshot | 64-138 |
| Selected links per snapshot | 1 at every snapshot |
| Distinct selected partner satellites | 13 |
| Median availability update | 0.334 ms |
| Median scoring | 0.00656 ms |
| Median matching | 0.0280 ms |
| Median total scheduling update | 0.369 ms |
| Maximum total scheduling update | 0.678 ms |
| Median allocated bytes per scheduling update | 160 bytes |

Measurements follow warm-up and exclude constellation generation, analytic
propagation, validation and CSV writing. They include the returned selection
vector allocation; workspace arrays are reused. Timings are observations on
this machine, not universal performance thresholds.

### Full-engine 4,000-spacecraft run

| Measurement | Observed |
|---|---|
| Final simulated time | 20.0 s |
| Wall time inside measured `run_simulation` | 27.37 s |
| Cumulative allocations inside that call | 3.475 GB |
| Accepted / rejected steps | 7 / 0 |
| Scheduler history samples | 8 (initialization + accepted steps) |
| Selected links | 1 at all 8 samples |
| Integrated target laser semimajor-axis benefit | +0.5359456026 m |
| Final spacecraft count / finite-state checks | 4,000 / PASS |

The full-engine measurement includes setup, state isolation and first-use
compilation in that process. Allocation is cumulative, **not peak resident
memory**. The relatively large engine allocation remains a limitation; this
work does not establish that long-duration 4,000-spacecraft propagation is cheap.

### Controlled matcher-only before/after measurements

For an N-satellite star, assign available edge `(1,j)` score `Float64(j)`, zero
penalty, warm up `select_interlinks!`, then measure one call with `@timed`.
All measurements selected the single highest-scoring edge.

| N | Original time | Original allocated bytes | Updated time | Updated allocated bytes |
|---|---:|---:|---:|---:|
| 128 | 13.26 ms | 24,244,072 | 0.00245 ms | 160 |
| 512 | 1,576.96 ms | 1,338,152,160 | 0.00715 ms | 160 |
| 4,000 | Not run | Not measured | 0.0549 ms | 160 |

These isolate matching, not availability or propagation. The original matcher
was not run at 4,000 nodes after its 512-node allocation reached 1.34 GB.
No before/after end-to-end speedup is claimed.

## Regression verification

**407/407 assertions passed**, no failures or skips:

- 221 array/matching assertions: signed transitions, eligibility/activity/battery
  changes, bulk and single-key APIs, workspace reuse, registration invalidation,
  60 small graphs compared with exhaustive enumeration, disconnected cycles,
  tie determinism independent of registration order, and 4,000-node star,
  chain and disjoint-edge graphs. Warmed bulk refresh and single-link tree
  matching are each checked against a 1,024-byte allocation ceiling.
- 43 Walker assertions: seeded reproducibility, satellite/candidate counts,
  circular geometry, plane/slot spacing, one-period analytic propagation,
  schedule output shape and invalid-input rejection.
- All 143 original assertions: terminal restrictions, physical scores and forces,
  accepted/rejected-step callback behavior, continuous-force integration,
  non-interlink behavior and the one-hour step-cap comparison.

The original 10 s versus 5 s step-cap comparison still gives maximum differences
of 0.146279331 m in integrated laser semimajor-axis change, 0.000138100594 m/s
in integrated laser velocity change, and 0.713948850 m in position, within its
existing limits. The persisted CSV files were also loaded back and their row
counts and 4,000-spacecraft result checked.

## Saved outputs

All experiment CSVs are in
[III_output/walker_interlinks/seed_20261003](III_output/walker_interlinks/seed_20261003/):

- [constellation.csv](III_output/walker_interlinks/seed_20261003/constellation.csv):
  complete generator settings.
- [snapshot_metrics.csv](III_output/walker_interlinks/seed_20261003/snapshot_metrics.csv):
  61 rows of component timings, allocations and link counts.
- [snapshot_schedule.csv](III_output/walker_interlinks/seed_20261003/snapshot_schedule.csv):
  selected endpoint pairs at snapshot times.
- [engine_summary.csv](III_output/walker_interlinks/seed_20261003/engine_summary.csv):
  full-engine timing and validation measurements.
- [accepted_step_schedule.csv](III_output/walker_interlinks/seed_20261003/accepted_step_schedule.csv):
  selected endpoint pairs from actual engine callbacks.
- [accepted_step_counts.csv](III_output/walker_interlinks/seed_20261003/accepted_step_counts.csv):
  all callback times and counts, including any empty selections in future runs.

## Scope and remaining limitations

- Forest components now have linear-time exact matching; arbitrary cyclic
  components still have exponential worst-case complexity. A polynomial-time
  weighted blossom solver would be needed for large general cyclic graphs.
  There is no silent greedy fallback.
- This tests 4,000 satellites and 3,999 target-star candidate links, **not**
  all 7,998,000 possible undirected satellite pairs or a constellation-wide
  objective. Satellite count alone does not describe graph difficulty.
- Availability still uses the existing range/resource/custom-eligibility
  model. No new Earth occultation, pointing, thermal or communications model
  was introduced. At this altitude, the experiment's 2,000 km range is below
  the equal-altitude Earth-horizon chord, but this is not a general LOS model.
- Scheduling remains at accepted step endpoints; crossings are not localized
  as continuous events. History storage grows with accepted steps.
- The snapshot benchmark covers an hour; the full-engine stress test covers
  only 20 seconds. Longer durations, dense graphs and more terminals need
  separate validation.
