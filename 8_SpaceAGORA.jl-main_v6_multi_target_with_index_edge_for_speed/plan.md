# V6 interlink performance plan

Date: 2026-10-05. Scope: this v6 directory only; legacy, v5 and v5 copy remain unchanged.

The completed performance plan below is historical. The 2026-10-13 graph
refactor supersedes its dictionary-API and state-mirror requirements:
`model.graph` owns vertices, edges, parameters and state arrays;
`edge_index_lookup` provides keyed access without connection wrappers.
Sorted edge order is now temporary scheduler data, and availability transitions
are computed only in diagnostics from before/after snapshots. See the
[current API](scalability-validation.md#current-graph-api-2026-10-13).

## Request and acceptance criteria

- **REQ-001:** Reduce lookup, allocation and force-evaluation overhead while retaining the dictionary lookup API and all three v6 policies.
- **REQ-002:** Cache terminal topology, preserve exact weighted matching, terminal exclusivity, continuous forces and accepted-step scheduling.
- **REQ-003:** Run the same PDF Case 1 configuration under legacy, v5 and v6; save numerical results and first-run/warmed timings.
- **REQ-004:** Generate a seeded randomized Walker-Delta constellation with exactly 4,000 satellites and save valid schedules at 0:60:3600 seconds. Do not run the copy's separate 20-second full-engine Walker simulation.

## Implementation steps

### P1: Indexed storage and forces

- [x] Preserve dictionary connection objects while adding indexed keys, references, scores and incident-edge lists.
- [x] Refresh Boolean availability and signed Int8 transitions in place; keep the single-key API.
- [x] Use fixed-size position differences and visit incident links during force evaluation.

### P2: Cached exact matching

- [x] Cache a 2-by-E terminal endpoint matrix, adjacency and traversal buffers; invalidate the cache on registration.
- [x] Benchmark the existing implementation before adding a tree fast path.
- [x] Use iterative DFS and exact BigInt tree dynamic programming; retain NetworkX Blossom for cyclic components.
- [x] Verify subnormal/near-tied scores, cache reuse, topology changes, mixed components, 4,000-node forests and existing physics regressions.

The baseline 4,000-terminal chain took 6.835 seconds per matching call, justifying the optional tree fast path. Tree and cyclic solvers use the same exact integer encoding of Float64 scores. Equal-optimum link choices can differ from Blossom-only results, but remain deterministic with fixed topology and backend.

### P3: Comparison and Walker experiment

- [x] Run PDF Case 1 for 10 initial target periods with J2, 21 spacecraft, 10-second maximum steps and 1e-12 tolerances in all three projects.
- [x] Repeat five warmed passes per process; separate first-use compilation and compare all SMA traces and schedules.
- [x] Run 61 analytic Walker snapshots at 60-second intervals with seed 20261003, 4,000 satellites, 80 planes, 550 km altitude and 53-degree inclination.
- [x] Save constellation metadata, per-satellite initial elements, schedules, timing/allocation metrics and validation results.
- [x] Summarize changes, results and limitations in scalability-validation.md.

## Constraints and validation

Topology changes must use `register_allowed_link!`; direct edits to graph topology arrays
or the key-to-index dictionary are unsupported. Registration updates topology eagerly;
the scheduler refreshes its workspace after growth, preserving existing scores and
initializing appended scores to zero. Raw scores belong to `policy.matching.score`,
not the graph. Scoring and selection must use the same policy; switching graphs or
discarding the workspace resets scores. Graph availability and active flags remain
separate from scheduler-owned scores.
The bulk availability API returns a fresh Boolean snapshot. Diagnostics compare it
with pre-update availability to count transitions; no transition vector is stored
on the graph. These counts are not a substitute for geometry or rescoring.

No dense satellite-by-satellite matrix, snapshot callback, new package, orbital-model change or legacy/v5 source edit is required. The Walker benchmark uses 3,999 target-to-helper candidates under fixed-target `:gve_sma`; it is not an all-pairs scalability test or force-coupled propagation. First-index policy and multi-terminal graphs are tested separately.

## Requirement mapping

| Requirement | Implementation | Evidence |
|---|---|---|
| REQ-001 | Interlink model, indexed policy loops, incident-edge forces | test/runtests.jl; test/benchmark_interlinks.jl |
| REQ-002 | Cached exact hybrid matcher | Exact/cached matching tests and existing simulation regressions |
| REQ-003 | II_examples/compare_pdf_case1.jl | III_output/pdf_case1_comparison |
| REQ-004 | II_examples/walker_interlinks.jl | III_output/walker_interlinks/seed_20261003 |