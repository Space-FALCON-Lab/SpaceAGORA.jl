# EDG adapter boundary preparation

This is a design proposal for a future external energy-depletion guidance (EDG)
adapter under the [optional adapter standard](adapter_contract.md). It creates no
public facade, package, contract version or new supported execution mode. Scientific
qualification and integration of the experimental controller remain prerequisites;
this document does not establish manuscript equivalence or authorize extraction.

The inventory below uses core commit
[`f9eb5f13`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/f9eb5f13cf2f19f4327abb358d22513960bef23f).
The experimental comparison uses draft
[PR251](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/251) at
[`6888f007`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/6888f007db34c14e9431bddd968442873e80cee8),
whose base is `Gen-EDG` at `8d86d5edf2cfd459d9467047b0ddd29b90328a02`, not this
core revision. These are source observations, not an accepted combined tree.
Refresh this comparison against the scientifically accepted revision before
implementing the interface.

## Existing ownership to preserve

Paths in this table are relative to core. The existing `EDGServices` module is
private: its name does not make it the versioned facade required by the standard.

| Responsibility | Current owner | Boundary decision |
|---|---|---|
| Public EDG configuration, state and guidance model identities | `src/gnc/guidance/target_energy_bracketing.jl` | Retain the existing exported type identities and construction path in core; do not define duplicate types in a companion. |
| Prediction, heat and structural constraints, targeting and command decisions | Six kernels included by `src/gnc/guidance/aerobraking/typed_edg/algorithms.jl` | Candidate domain implementation for later extraction, after scientific acceptance. |
| Physical-state projections, atmosphere and frame queries | `src/gnc/guidance/aerobraking/typed_edg/services.jl` | Keep implementation in core; expose only reviewed value projections and query operations needed by the accepted algorithm. |
| Model callbacks and panel actuation | `src/gnc/control/targeting_control.jl` and the guidance model wrapper | Core invokes the selected algorithm and applies its validated command to configured panel links. |
| Fresh-run initialization and passage resets | `src/gnc/guidance/aerobraking/energy_depletion_lifecycle.jl`, `src/gnc/control/energy_depletion_lifecycle.jl` | Core schedules lifecycle events; one owner updates each shared EDG runtime. |
| Integrator, physical heat integrals, event roots, atmosphere mask, solver tolerances and output saving | `src/simulation/engine/` and `src/simulation/callbacks/` | Remain core-owned. The adapter must not mutate their storage or install replacement phase handlers. |
| Older E-EDG/T-EDG strategy and bridge paths | `src/gnc/guidance/aerobraking/e_edg/`, `t_edg/` and control bridge code | Outside the initial extraction scope; inventory their callers before changing any shared helper. |

The present algorithm still consumes `ODEParams`, `u`, spacecraft objects and
private core modules. Re-exporting that whole module would expose mutable engine
storage. The future facade must replace those accesses deliberately, with parity
tests, rather than treating the current signatures as a stable external API.

## Proposed supported interface

Use one guidance-domain facade with the standard's version check, compatibility
error and explicit implemented/consumed inventories. Its final name, symbols and
version are to be selected in the implementation review. The operations below
describe the proposed contract, not callable methods available today.

| Operation | Inputs | Result and permitted effects |
|---|---|---|
| Preflight | EDG configuration, spacecraft/link mapping, callback periods and execution capabilities | Validate supported combinations before output creation or stepping; no scientific feasibility claim. |
| Initialize | Run-owned EDG runtime and initial observation | Reset plan/history according to the reconciled lifecycle; establish any required passage baselines without resetting physical integrals. |
| Passage transition | Spacecraft index, event time, boundary kind, direction and copied observation | Update only that spacecraft's plan/history; deliver once to the runtime shared by guidance and control. |
| Guidance update | Immutable observation, configuration, run-owned EDG runtime and scoped prediction services | Select mode and update bracket/target state; return a diagnostic snapshot. No link rotation or trajectory-state write. |
| Control update | The same ownership inputs at the configured control tick | Return an apply flag, the core-owned panel-AoA command, switch/cache status and diagnostic snapshot; core applies the command. |

The control result should preserve the distinction already present in
`EDGControlDecision`: whether to apply a command, the command itself and what the
switch planner did. A successful computation or cache hit is not evidence that a
trajectory meets its constraints. Reuse the existing core command identity where
appropriate; do not return a mutable view of runtime state.

An observation must identify the 1-based spacecraft and controlled links, elapsed
simulation time in seconds, inertial position in metres, inertial velocity in
metres/second, mass in kilograms and copied per-panel cumulative heat in J/cm².
Environment observations must label geodetic altitude in metres, density in
kg/m³, temperature in kelvin, airspeed in metres/second, dimensionless molecular
speed ratio and dynamic pressure in pascals. AoA is in radians; heat rate is
W/cm²; specific energy is J/kg. Missing values and invalid samples need explicit
status, preserving the accepted algorithm's existing handling rather than
silently inventing replacements.

Core must distinguish simulation time from ephemeris time. Existing live and
trajectory queries use the configured epoch plus elapsed time for frames and the
elapsed time for atmosphere queries. Local atmosphere wind is east/north/up;
airspeed is spacecraft velocity relative to the rotating planet minus atmosphere
velocity. `l_pi` maps inertial components to planet-fixed components. Preserve
these conventions and the rotation transport term. The altitude-only prediction
helper currently queries zero latitude/longitude; it is not interchangeable with
the full trajectory query.

Prediction services should take explicit state/time values and a run-scoped,
opaque environment context. Inventory frame/atmosphere queries, configured gravity
evaluation, per-link aerodynamic evaluation and vehicle geometry projections.
Retain core implementations of shared physics; do not copy them into the adapter.
The facade must specify query side effects, atmosphere caches and RNG ownership:
a predictor must not gain direct access to `shared_buffers`, mutate geometry, or
silently reseed a live atmosphere. Separate prediction instances or replayable
samples require a reviewed equivalence decision if they change existing query
semantics. Call-scoped borrowed service handles must not escape their run.

## State ownership and lifecycle

Core retains the public configuration/model/state types and creates the run's
configuration. The EDG implementation may mutate only its run-owned algorithm
state through the agreed contract. Paired guidance and control must share that
same state within a run; separate runs and supported concurrent samples must not
share it. Current `run_simulation` isolation deep-copies configuration while
preserving that within-run relationship. Explicit `isolate_state=false` means
caller ownership and fresh initialization, not implicit resume.

Current preflight checks require state-vector lengths to match spacecraft count,
at most one EDG guidance and one EDG control model, matching paired configuration
and shared state, and a guidance partner for targeting control. Maximum-depletion
control can operate alone and then owns initialization/reset. Preserve these
guards unless a separately reviewed feature deliberately changes them. Record
any additional adapter scratch arrays and RNGs as per-run resources; define
cleanup if a later implementation allocates external resources. There is no
need for a global provider registry or a dynamics-owner stepping loop.

The engine owns event scheduling and normal guidance/control periods. It supplies
the correct spacecraft index and calls the adapter at those times; the adapter
does not advance the physical trajectory. Core alone rotates the validated
non-root panel links. EDG's current force/torque effector returns zero; the normal
aerodynamic dynamics observe the resulting geometry. Preserve callback ordering,
held-command activation, event save behavior and output sampling explicitly in
the accepted integration. No period, threshold or solver setting changes as a
side effect of packaging.

### Main/experimental differences that must be reconciled

| Topic | Inspected core | Experimental PR251 observation | Required decision before facade implementation |
|---|---|---|---|
| Passage plan lifetime | The engine's atmosphere transition invalidates cached plans, retaining counters/last telemetry; duplicate unchanged masks do not reset again. | Branch bookkeeping adds per-panel geodetic entry/exit heat snapshots and retains its branch's scheduling policy. | Specify exactly which boundary invalidates a plan and which captures/finalizes heat, including coincident events, duplicate delivery, start-inside, poles and multiple passages. Do not equate the spherical phase boundary with the geodetic EDG gate. |
| Heat quantity | The typed decision path reads cumulative controlled-link heat. | Planning uses the maximum of individual panel increments for the current passage; physical heat and existing telemetry remain cumulative. | Qualify the passage accounting and its transfer to current main. Never subtract a maximum baseline from a maximum current load, or reset physical heat integrals. |
| Checkpoints | Preflight rejects both checkpoint writing and resume because EDG history is absent. | Initialization explicitly rejects nonzero-time resume without passage baselines/history. | Keep current core's refusal until a separately versioned complete history format and continuation tests exist; this proposal grants no checkpoint support. |
| Scientific policy | Current typed kernels and lifecycle are the executable baseline. | Applied-profile, terminal/event, passage-constraint and sensitivity acceptance remain open. | Select one reviewed scientific baseline, integrate it with current main and retain discriminators before moving code. Do not choose equations or cadence in this document. |

## Error behavior

An incompatible contract must fail before adapter methods or optional libraries
become available, and be checked again during runtime initialization. Missing
providers and unsupported capabilities must produce a clear setup error, not
activate a different algorithm. Invalid configuration, ownership or panel mapping
must fail before actuation. Preserve current checkpoint refusal.

Prediction failure, an unreachable target, a sampled constraint exceedance and
absence of a valid command are distinct outcomes. Preserve the scientifically
accepted mode fallback and command law, accompanied by structured diagnostics;
do not introduce a new generic low-drag override or classify a sampled failure as
proof of global infeasibility. Programmer errors must remain errors. Existing
legacy exception fallback behavior must be inventoried before any change in
propagation policy. A rejected compatibility load requires a fresh process.

## Acceptance sequence and evidence

| Stage | Required evidence | Status of this preparation |
|---|---|---|
| Scientific baseline | Exact accepted source, cases, settings and manuscript interpretation; applied commands, terminal/event certificates, passage constraints and sensitivity evidence. Retain historical results with their original labels. | Pending separate scientific review; no new acceptance here. |
| Current-main reconciliation | Explicit resolution of the differences above; source review and unchanged scientific tolerances; existing ownership, lifecycle, wind-relative and strategy-parity regressions plus the accepted branch's heat/predictor discriminators. | Not performed by this document. |
| Public boundary | Explicit bindings, exact version refusal, no private engine/storage access; state/time/unit fixtures and unsupported-mode errors. Extend `test/unit/adapter_contract_tests.jl` to the new facade only when implemented. | Proposed. Existing shared tests already cover HYPR and external propagation. |
| Lifecycle and isolation | Fresh/reused configurations, shared paired state, control-only case, two spacecraft, independent concurrent runs, duplicate/coincident transitions, start-inside/at-boundary, two passages and checkpoint refusal. Verify physical heat and noncontrolled geometry stay core-owned. | Use current lifecycle tests as a baseline; add reconciled passage cases. |
| External installation | Declared dependencies and exact qualified source pair; fresh core-first/adapter-first loading, repeated initialization, early incompatibility refusal, core-only isolation and one end-to-end configuration/run/result test. Require native validation for any supported native atmosphere path. | Future adapter work. |
| Extraction | Publish reviewed core facade, qualify adapter against it, then pin that adapter in mandatory current-core CI. Preserve type identities, source history, licenses and fixtures; remove the in-tree implementation only after both repositories pass. | Blocked on preceding stages; no extraction authorized by this design. |

Tests of module shape establish interface conformance, not scientific validity.
The initial scope excludes moving the legacy bridge paths, redesigning EDG
equations, changing controller settings, enabling checkpoints or claiming new
concurrency support. Those require their own evidence and review.
