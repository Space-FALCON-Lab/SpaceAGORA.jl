# GNC Aerobraking Boundary Contract

## Canonical ownership
- `src/gnc/guidance/aerobraking/*` owns aerobraking strategy planning logic (`E-EDG`, `T-EDG`).
- `src/gnc/control/aerobraking/*` owns algorithm-agnostic tracking/execution only.
- `src/mission/operations/aerobraking_policy/*` owns strategy selection policy interfaces.

## Required boundaries
- Guidance files must not include control source files directly.
- Control files must not define targeting/switch-planning strategy solvers.
- Strategy selection must be typed and mission-owned (`AerobrakingStrategyKind`, selector contract).
- DRL integration point is a stub contract only (`DRLPolicyAdapterStub`), with explicit `Not implemented` behavior.
- Propulsive maneuver coupling must cross the typed command boundary (`PropulsiveManeuverCommand`, `AerobrakingControlCommand`) rather than `DynamicEffectors` or direct thruster-state mutation from guidance.

## Canonical interfaces
- `AerobrakingGuidanceInput`
- `AerobrakingGuidanceOutput`
- `PropulsiveManeuverCommand`
- `AerobrakingControlCommand`
- `AbstractAerobrakingStrategy`
- `AerobrakingPolicyConfig`
- `AbstractAerobrakingPolicySelector`
- `AerobrakingStrategyKind` (`E_EDG`, `T_EDG`)


## Typed EDG prediction and panel commands

The typed energy-depletion route uses one calculation owner,
SimulationModel.EDGAlgorithms, under
src/gnc/guidance/aerobraking/typed_edg/. It owns reachable-energy brackets,
mode decisions, heat-load windows, targeting switches, heat/structural limits,
and the constrained angle decision. The private EDGServices module supplies
SpaceAGORA state, atmosphere and frame projections. ODEParams remains an
internal dependency; this change does not establish a standalone EDG package.

Guidance and control model types remain in GuidanceHooks and ControlHooks.
Their public constructors, fields and hook generics retain their owners.
The guidance hook forwards at its existing scheduled times. The control hook
requests control_decision! at its existing control times, then applies the
returned AerobrakingControlCommand to its panel effector exactly once.
Control-first scheduling and the accepted pass lifecycle are unchanged.

The decision takes config, shared state, explicit actuator/measurement links,
the current state, ODEParams, callback seconds and the one-based spacecraft
index. Prediction still uses config.controlled_panel_links, which may differ
from the control effector's links. The result contains an angle in radians,
a separate apply flag, a switch_action (:cached, :outside_pass, :solved,
:not_required or :invalid_index), and already-computed telemetry.
The :solved status records execution, not scientific acceptance.
An invalid index produces no command application or environment query.

Numerical kernels do not move panels or propagated heat loads. The decision
updates the existing EDG cache and telemetry before panel application, in the
same order as the former control hook. Existing exceptions still propagate;
a decision can update mode before an atmosphere error, and completed telemetry
can precede an actuator error. There is no rollback or new error suppression.

Targeting samples current latitude/longitude and uses spacecraft velocity
minus local ENU wind transformed to planet-fixed coordinates. The heat-load
sampler retains its altitude-only, zero-latitude/longitude query and ignores
wind. Frame queries use et_start + callback_seconds. These routes, lazy early
returns, query order/count, absolute switch times, inclusive heat-load
endpoints, certification reasons and Inf sentinels are preserved.
The recomputation interval remains inactive.

Old qualified control numerical names are imported aliases to the single
owner; their signatures and keyword defaults travel with their definitions.
The model-taking private control methods are explicit forwarding methods.
Guidance's private numerical compatibility methods forward to the later-loaded
owner. Their former private generic identities are not promised across the
refactor. Duplicate state projections are aliases to EDGServices. Legacy
E-EDG/T-EDG wrappers and their bridges retain their existing routes.

The required GNC boundary gate also checks this typed owner for duplicate
definitions, reverse dependencies on control and panel mutation. The focused
ownership tests cover command/application parity, invalid-index laziness,
switch endpoints, input forms, panel-selection mismatch, query coordinates,
epoch and wind conventions, plus deliberately broken ownership, wind-sign and
time-offset controls. Numerical equivalence is assessed on a matched parent
and candidate with the retained bounded predictor fingerprint; this does not
replace full mission or scientific acceptance.
