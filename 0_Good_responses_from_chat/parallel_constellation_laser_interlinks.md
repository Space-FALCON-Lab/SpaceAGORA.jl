# Parallelism for Constellations with Laser Interlinks

## The Key Distinction

A constellation with laser interlinks is a **coupled system**: the force on one spacecraft can depend on another spacecraft's current position. The solver therefore needs to integrate the combined constellation state together with a full-state solver.

This does **not** mean that all computation must be serial. Parallel execution can still be used inside each right-hand-side (RHS) evaluation. A full-state solver and parallel RHS evaluation describe different things:

- **Full-state solver:** advances the coupled state of all spacecraft together.
- **Parallel RHS execution:** calculates parts of the derivative concurrently from one shared state snapshot.

Split solver modes such as `:gravity_backbone_split`, `:split_imex`, and `:multirate` are currently rejected when interlinks are enabled. This restriction does not prohibit the regular parallel RHS execution modes.

## What Happens During an RHS Evaluation

At one solver evaluation, the solver provides a single current state snapshot, `u(t)`, containing every spacecraft. Workers read from that snapshot and compute derivatives; no worker advances its spacecraft state in place.

Conceptually:

```text
Shared current state u(t)
    Worker 1 reads relevant spacecraft positions and computes assigned derivatives
    Worker 2 reads relevant spacecraft positions and computes assigned derivatives
    Worker 3 reads relevant spacecraft positions and computes assigned derivatives
All calculations finish
Solver uses the complete derivative to calculate its next stage or state
```

For a laser link between spacecraft `i` and `j`, the force direction is based on their positions in the same snapshot:

```text
separation = position[i] - position[j]
direction = separation / norm(separation)
```

The direction is recalculated for each RHS evaluation as the solver evaluates new stage states. It is not taken from a worker's partially advanced or stale spacecraft state.

## Example: Four Spacecraft

Suppose four spacecraft have two active links:

```text
SC1 <-> SC2
SC3 <-> SC4
```

If the interlink pass were parallelized by spacecraft, separate workers could calculate the per-spacecraft contributions:

| Worker | Reads | Calculates and writes |
| --- | --- | --- |
| 1 | SC1 and SC2 positions | Force on SC1; writes only SC1's derivative |
| 2 | SC2 and SC1 positions | Force on SC2; writes only SC2's derivative |
| 3 | SC3 and SC4 positions | Force on SC3; writes only SC3's derivative |
| 4 | SC4 and SC3 positions | Force on SC4; writes only SC4's derivative |

Each worker can read the complete shared state while writing to a distinct derivative slot. This avoids races. The solver waits for the derivative calculations before progressing.

If a spacecraft can have multiple active links, its worker can sum the forces from its neighbors and then write its derivative once. Do not parallelize by link if that makes multiple workers scatter writes into the same spacecraft derivative.

## RHS Execution Modes

The `plan.mode` field selects an RHS execution strategy; it is not the ODE solver mode.

- `:serial`: performs spacecraft and effector calculations in sequence on one worker.
- `:satellite_batch`: distributes spacecraft across workers; each spacecraft's effectors are calculated in sequence by its worker.
- `:per_satellite_effector_reduce`: distributes a spacecraft's effector calculations across workers, then sums their force and torque contributions.
- `:flat_constellation_effector_queue`: schedules satellite/effector work items through a flat constellation-scale queue and reduces per-worker partial results.

### Comparing the Modes: Three Spacecraft, Three Effectors

Suppose SC1, SC2, and SC3 each have gravity, solar-pressure, and drag effectors. That gives nine ordinary spacecraft/effector calculations for one RHS evaluation.

| Mode | Example distribution |
| --- | --- |
| `:serial` | One worker calculates all three effectors for SC1, then SC2, then SC3: nine calculations in sequence. |
| `:satellite_batch` | Up to three spacecraft jobs run concurrently. A worker assigned SC1 calculates its gravity, solar-pressure, and drag effects in sequence. |
| `:per_satellite_effector_reduce` | The effectors for a spacecraft, for example SC1's gravity, solar pressure, and drag, can be calculated concurrently; their force and torque results are then summed. The planner favors this inner-effector strategy when it is a better fit than dividing work primarily by spacecraft. |
| `:flat_constellation_effector_queue` | Work is exposed as items such as `(SC1, gravity)`, `(SC1, solar pressure)`, and `(SC2, drag)`. Workers take available items, and partial results are combined by spacecraft afterward. |

Nine queue items do **not** mean nine workers. With four available workers, at most four eligible items can be in flight at once. Flat queuing can balance uneven workloads, such as one spacecraft having a much slower effector, but it adds job-dispatch, synchronization, temporary-result, and reduction overhead. When the work is light or evenly balanced, `:satellite_batch` can be as fast or faster because each worker keeps its spacecraft's results locally. Performance depends on measured work and configuration; flat execution is not always quicker.

The planner can choose among these strategies based on settings and runtime conditions. `SPACEAGORA_RHS_EXECUTION_MODE=auto` requests automatic selection; `auto` is not itself a `plan.mode` value.

The flat queue is one parallel strategy, not a requirement for a constellation. Both it and satellite-batch execution work within a full-state solver.

These are RHS strategies for the ordinary dynamic effectors. In the current implementation, `_apply_interlink_rhs!` adds laser interlink forces in a separate serial pass; it is not automatically split into the flat queue's satellite/effector jobs.

## Scheduling Versus Force Evaluation

The scheduler chooses which candidate links are active at initialization and at its discrete callback points. Between scheduling updates, the active link topology is held fixed. However, for each RHS evaluation, the force direction and magnitude for those active links are recalculated using that evaluation's current spacecraft positions.

Thus the **active links** usually change less often than the **force values**.

## Current Implementation and Possible Parallelization

Currently, regular dynamic-effector work can use parallel RHS modes, while `_apply_interlink_rhs!` applies interlink forces in a separate serial pass. The interlink force code reads the current full state and adds the result to spacecraft derivatives.

A race-free parallel interlink pass could divide work by spacecraft: each worker reads the current shared state and active-link graph, sums forces for one spacecraft, and writes only that spacecraft's derivative. An adjacency list of active neighbors could avoid scanning every candidate link for every spacecraft. Parallelizing by link would require a reduction or other synchronization because multiple links can contribute to the same spacecraft.

Whether this improves runtime depends on constellation size, active-link count, and thread-dispatch overhead. The current serial implementation is simpler for small cases; a parallel version should be benchmarked and checked for trajectory parity before adoption.

## Extra Note:

so between the calback for multiple internal step for propogation, are each worker get the force direction from the satellite positions from the last internal step, or the updated global satellite positions within the current step? if it's the later, then within one step sokme satellite might propagate early and some might finish propogation later, and their might be mismatch?

11:43 PM
It’s the second idea, but there is no “some satellites propagate early” within a step. The solver does not update each spacecraft’s position as its worker finishes.

At each RHS evaluation, the solver passes one shared state snapshot 
𝑢
(
𝑡
)
u(t) to the calculation. Workers read positions from that same snapshot, compute their assigned derivatives, and write to separate slots in 
𝑢
˙
u
˙
 . Once those calculations are done, the solver uses the complete derivative to calculate its next internal stage or step state.

For example, when computing SC1’s laser-force direction, its calculation uses SC1 and SC2 positions from the same 
𝑢
(
𝑡
)
u(t), even if another worker is calculating SC2’s derivative. No worker advances the state itself, so there’s no mid-step position mismatch. The callback may change which links are active at its scheduling time; the force direction is then recomputed from the current snapshot at each RHS evaluation.

