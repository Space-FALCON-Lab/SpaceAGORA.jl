# Main differences: main_v5, main_v5 copy, and main_v6

This document compares the current source in:

- [main_v5](../8_SpaceAGORA.jl-main_v5/)
- [main_v5 copy](<../8_SpaceAGORA.jl-main_v5 copy/>)
- [main_v6](../8_SpaceAGORA.jl-main_v6/)

These are different implementations, not identical copies. The examples below
use hypothetical scores to explain their behavior.

## Comparison table

| Area | main_v5 | main_v5 copy | main_v6 |
|---|---|---|---|
| Scheduling policies | SMA and eccentricity | SMA only | SMA, eccentricity, and first-index SMA |
| Target selection | Specific target, or `target_idx=nothing` to sum contributions from both endpoints | One specific target only | Same as v5, plus per-pair lower-index targeting |
| Matching algorithm | NetworkX weighted Blossom for all eligible graphs | Tree dynamic programming; recursive subset matching for cyclic components | Tree dynamic programming; NetworkX Blossom for cyclic components |
| Matching arithmetic | Exact integer weights converted from Float64 scores | Float64 score sums | Exact integer weights for both tree and Blossom solvers |
| Storage | Primarily dictionary-based | Dictionary plus indexed edges, availability arrays, and cached matching workspace | Similar indexed structure, also including a score array |
| Force evaluation | Scans all candidate links for each spacecraft | Visits only that spacecraft's incident links | Visits only that spacecraft's incident links |
| Custom eligibility hook | No | Yes | No |
| Per-selected-link penalty | No | Yes | No |
| Matcher dependency | PythonCall + NetworkX | No PythonCall/NetworkX requirement for matching | PythonCall + NetworkX |
| Examples/benchmarks | Small GVE example, PDF Case 1, matching benchmark | Walker example and scalability tests | PDF comparisons, Walker snapshots/integration comparisons, interlink phase benchmarks |

## 1. Scheduling policies: how useful is a link?

A scheduling policy assigns a score to each candidate link.

- **SMA** means semimajor axis: a measure of orbital size. The SMA policy scores
  the instantaneous rate at which a laser force increases that size.
- **Eccentricity** describes how noncircular an orbit is. This policy scores the
  instantaneous eccentricity increase caused by the force; it is not a policy
  for reducing eccentricity or tracking a desired eccentricity.
- **First-index SMA** scores the SMA increase of the lower-index spacecraft in
  each pair.

For example, first-index SMA scores spacecraft 3 for link 3-8, but spacecraft 6
for link 6-9. It does not have one global target.

These scores are predicted instantaneous rates, not current orbital values or
guaranteed final orbit changes. Nonpositive scores are excluded from matching;
the copy can additionally subtract a penalty.

Sources: [v5 policies](../8_SpaceAGORA.jl-main_v5/I_src/1.4_satellite_interlink/scheduling_policies.jl),
[copy policy](<../8_SpaceAGORA.jl-main_v5 copy/I_src/1.4_satellite_interlink/scheduling_policies.jl>),
[v6 policies](../8_SpaceAGORA.jl-main_v6/I_src/1.4_satellite_interlink/scheduling_policies.jl).

## 2. Target selection: whose benefit is counted?

With ordinary SMA and `target_idx=3`, only spacecraft 3's predicted SMA rate
counts. Links not involving spacecraft 3 get zero score.

In v5 and v6, `target_idx=nothing` instead adds the rates of both endpoints.
The copy requires one designated target.

V6's first-index SMA ignores `target_idx` and chooses the lower spacecraft
vector index separately for every link. This index is not altitude,
registration order, or the spacecraft's separate `id` field.

Suppose link 1-3 gives spacecraft 1 a rate of +2 m/s and spacecraft 3 a rate
of -1 m/s:

| Configuration | Score for link 1-3 |
|---|---:|
| Ordinary SMA, `target_idx=3` | -1 m/s: excluded |
| Ordinary SMA, `target_idx=nothing` | +1 m/s: both rates added |
| First-index SMA | +2 m/s: only spacecraft 1's rate counted |

First-index SMA does not require both endpoints to benefit.

## 3. Matching algorithm: which links can operate together?

A matching selects links without reusing a terminal. With one terminal per
spacecraft, a spacecraft cannot participate in two selected links. With
multiple terminals, the restriction is per terminal, not per spacecraft.

Consider this chain, where the numbers are scores:

```text
A --6-- B --11-- C --6-- D
```

Selecting B-C gives 11. Selecting A-B and C-D gives 12, so the two-link
selection wins. Matching is not simply choosing the largest individual score.

The versions solve this problem differently:

- **V5:** uses NetworkX's weighted Blossom solver, which handles general graphs.
- **Copy:** splits the graph into disconnected components. Components without
  loops use tree dynamic programming, which combines smaller solutions from
  the ends inward. Components with loops use recursive subset matching with
  remembered results.
- **V6:** uses tree dynamic programming for components without loops and
  Blossom for components with loops.

A chain is a tree. A triangle A-B-C-A contains a cycle.
The classification uses currently available, positive-weight links.

The copy's cyclic subset search can require exponentially many combinations,
so a large cyclic network can become very expensive. V6 avoids that approach
by using polynomial-time Blossom for cyclic components. Tree matching is
linear in component size in arithmetic operations; exact large-integer
operations in v6 have their own additional cost.

Sources: [v5 matcher](../8_SpaceAGORA.jl-main_v5/I_src/1.4_satellite_interlink/scheduling_models.jl),
[copy matcher](<../8_SpaceAGORA.jl-main_v5 copy/I_src/1.4_satellite_interlink/scheduling_policies.jl>),
[v6 matcher](../8_SpaceAGORA.jl-main_v6/I_src/1.4_satellite_interlink/scheduling_models.jl).

## 4. Matching arithmetic: how are totals compared?

The copy adds and compares `Float64` scores. These are ordinary scientific
floating-point numbers with limited precision. Very close totals can be
affected by rounding.

V5 and v6 convert the stored Float64 scores into exact fractions, then scale
them into `BigInt` integer weights:

```text
Original scores:       0.5, 0.25, 0.75
Multiply all by 4:       2,    1,    3
```

Multiplying every score by the same positive number preserves which total is
largest. `BigInt` permits exact addition even when the integers become large.
The copy can round when adding scores; v5 and v6 do not round those integer
matching sums.

This makes matching comparisons exact for the already-computed scores. It
does not make the physics or numerical integration exact, or restore precision
lost during score calculation. Equal-scoring solutions can also have different
selected links across solvers.

## 5. Storage: how is link information organized?

V5 primarily uses a dictionary: a link key retrieves a connection object.

```text
A-B -> connection with availability, score, and active state
B-C -> connection with availability, score, and active state
C-D -> connection with availability, score, and active state
```

The copy and v6 retain that dictionary and also number the links:

| Edge index | Link | Available |
|---|---|---|
| 1 | A-B | true |
| 2 | B-C | true |
| 3 | C-D | false |

They store which edges touch each spacecraft:

```text
A -> [1]
B -> [1, 2]
C -> [2, 3]
D -> [3]
```

Both also keep availability arrays, signed availability transitions
(`+1` entered, `-1` left, `0` unchanged), and reusable matching workspace.
V6 additionally keeps a dedicated score array.

This is like keeping both a contact directory and a numbered list optimized
for repeated processing. It needs additional bookkeeping but avoids repeated
lookups and reconstruction. The dictionary and indexed connection list refer
to the same connection objects; they are not independent copies of the physics.

These structures store registered edges, not a full spacecraft-by-spacecraft
matrix. In the indexed versions, topology changes must use
`register_candidate!` so the lists and matching cache stay consistent.

Sources: [v5 storage](../8_SpaceAGORA.jl-main_v5/I_src/1.4_satellite_interlink/inter_link_models.jl),
[copy storage](<../8_SpaceAGORA.jl-main_v5 copy/I_src/1.4_satellite_interlink/inter_link_models.jl>),
[v6 storage](../8_SpaceAGORA.jl-main_v6/I_src/1.4_satellite_interlink/inter_link_models.jl).

## 6. Force evaluation: which links are inspected?

To compute spacecraft A's laser force:

- **V5** scans every candidate link, looking for active links involving A.
- **Copy and v6** inspect only A's stored incident-edge list.

In the storage example above, A needs only edge 1. There is no need to inspect
edges 2 and 3.

Across N spacecraft and E candidate links, this reduces edge inspections from
O(N*E) to O(N+E). This is an efficiency change, not a different laser force law.
It does not imply the whole simulation has the same speedup: integration and
other dynamics still require work.

## 7. Custom eligibility hook: an extra user-defined rule

Only the copy lets you supply an extra function that returns true or false
for a candidate link.

For example:

> Reject this link if Earth blocks the line of sight, even when it is within range.

The built-in checks run first, including range, terminal validity, battery,
temperature, activity, and forbidden pairs. If they pass, the optional custom
function can reject the link.

The hook does not automatically implement Earth blockage; the user must supply
that logic. It cannot override a failed built-in check. Its default is `nothing`,
meaning no extra rule.

Source: [copy availability logic](<../8_SpaceAGORA.jl-main_v5 copy/I_src/1.4_satellite_interlink/inter_link_models.jl>).

## 8. Per-selected-link penalty: a fixed cost for using each link

Only the copy provides `active_link_penalty`. Its matching objective is:

```text
Total score - penalty * number of selected links
```

Using the earlier chain, with a penalty of 2:

| Selection | Raw total | After penalty |
|---|---:|---:|
| B-C | 11 | 11 - 2 = 9 |
| A-B and C-D | 12 | 12 - 4 = 8 |

Without a penalty, two links win. With this penalty, one link wins.

The default penalty is zero. Links scoring at or below the penalty are excluded.
The penalty has the same units as the score.

This is not a switching penalty: it does not specifically discourage changing
partners. It also does not itself model battery consumption or alter the force.

Source: [copy matching weights](<../8_SpaceAGORA.jl-main_v5 copy/I_src/1.4_satellite_interlink/scheduling_policies.jl>).

## 9. Matcher dependency: external software needed for matching

- **PythonCall** allows Julia to call Python.
- **NetworkX** is a Python graph library providing the Blossom matcher.

V5 uses both for nonempty eligible matching problems.
V6 solves tree components in Julia and uses NetworkX for cyclic components.
It still declares PythonCall as a project dependency.
The copy's matcher is implemented in Julia without PythonCall or NetworkX.

Analogy: v5 sends matching work to an external specialist; the copy handles
it internally; v6 handles trees internally and sends cyclic work to the specialist.

V5 and v6 configure NetworkX through CondaPkg. First-use provisioning can
require network access. The copy's lack of a Python matcher dependency does
not mean that its whole simulation has no other dependencies.

## 10. Examples and benchmarks: what can be demonstrated or measured?

- **V5:** a small GVE interlink example, the PDF Case 1 comparison, and a
  matching benchmark.
- **Copy:** a Walker constellation example and dedicated scalability and
  Walker tests.
- **V6:** PDF comparisons, Walker snapshot and force-coupled integration
  comparisons, and a benchmark measuring availability, scoring, matching,
  and force evaluation separately.

A snapshot scheduling benchmark measures scheduling at supplied positions
and velocities. It is not the same as integrating laser-driven orbital motion.
A matcher-only speedup is likewise not a whole-simulation speedup.

For recorded results and limitations, see the
[copy scalability report](<../8_SpaceAGORA.jl-main_v5 copy/scalability-validation.md>)
and [v6 scalability report](../8_SpaceAGORA.jl-main_v6/scalability-validation.md).
This comparison document does not represent a new test or benchmark run.

## Practical summary

- **V5:** simpler dictionary-based baseline with general-purpose Blossom matching.
- **Copy:** efficient indexed storage and tree matching, plus eligibility and
  penalty options, but fewer policies and an expensive cyclic subset solver.
- **V6:** indexed storage and tree optimization combined with Blossom for cycles,
  exact matching arithmetic, and first-index SMA.

V6 is the strongest starting point for general link networks, but it is not a
strict superset of the copy: the eligibility hook and link penalty are absent.
