# V5 Laser-Link Reconstruction Plan


## 1. SpacecraftModel

Add the following fields to `SpacecraftModel`:

| Field | Meaning |
| --- | --- |
| `n_terminal` | Number of physical connection terminals on the satellite |
| `health_index` | Satellite health on a scale from 0 to 100 |
| `tempurature_index` |Overall tempurature-wise health for satellite subsystems on a scale from 0 to 100 |

Each terminal can participate in at most one active connection at a time.
A satellite with several terminals may therefore have several simultaneous
connections, subject to other eligibility constraints.

Note: The existing spacecraft `Link` type represents a mechanical body. Keep that
meaning separate from intersatellite connections. The inter-satellite links is therefore referred to as "interlink" instead of "link"

## 2. InterLink Model and Terminal Graph

Introduce a mutable `InterLinkModel` containing a link type, such as laser or radio,
and a typed dictionary of candidate connections. Radio is a possible type in
the design, not a request to implement radio physics now.

Use terminal endpoints rather than satellite-only keys:

```julia
linkgraph[((satellite_i, terminal_m), (satellite_j, terminal_n))]
```

Here, `(satellite_i, terminal_m)` identifies terminal `m` on satellite `i`.
Each dictionary value contains the connection's parameters and mutable state.
Concrete type and function spellings remain implementation choices.

### Link Parameters

Keep physical configuration separate from scheduling state:

| Parameter | Meaning |
| --- | --- |
| `P` | Laser power, in watts |
| `B` | Magnification factor |
| `range` | Permitted connection range, in meters |

### Link State

Use a mutable state record with:

| State | Type | Initial value |
| --- | --- | --- |
| `available` | `Bool` | `false` |
| `active` | `Bool` | `false` |
| `score` | `Float64` | `0.0` |

Using `Bool` replaces the originally proposed integer 0/1 flags.

The dictionary may start empty during construction, but candidate terminal
pairs and their parameters must be registered before availability updates can
operate on them. Registering a candidate does not make it available or active.

## 4. Availability Checks

An availability function evaluates a candidate connection using the current
satellite states and configuration. It checks:

1. The endpoints are within the permitted range.
2. Both satellites satisfy the configured health threshold.
3. The connection is permitted by restrictions such as forbidden satellite pairs.
4. Any additional explicitly configured eligibility constraints are satisfied.

Return a Boolean and update the candidate's `available` field. Availability
means eligible for selection; it does not mean selected.

Terminal occupancy is a constraint on the selected set. Do not make a candidate
ineligible merely because a competing connection occupied its terminal in the
previous schedule, since the optimizer must be able to replace that connection.

## 5. Scheduling Policy and Scores

Introduce a separate simulation-level `SchedulingPolicyModel`, alongside the
existing guidance, navigation, and control configuration concepts. The requested
minimal model stores the policy selection rather than duplicating that selection
inside every physical link.

A scoring function evaluates each available candidate from a shared current
state snapshot and the selected policy. For example, a policy may score a desired
change in an orbital element using Gauss variational equations.

The discussed starting approach is an instantaneous score, not a prediction
multiplied by the unknown next adaptive solver-step duration. The scheduler does
not need to know that future duration to rank instantaneous benefits.

Scores must represent a consistent objective with a consistent sign and units.
For a constellation-wide objective, define how the effects on the participating
satellites contribute. A deorbit policy, for example, should reward the desired
decrease rather than blindly maximize a signed semimajor-axis derivative.

When the score represents physical orbital benefit, use the actual acceleration
and appropriate mass and power dependence. An unweighted sensitivity to a unit
force direction is not automatically comparable across unequal spacecraft or
links. The exact score formula and constellation aggregation remain to be defined.

## 6. Optimal Selection and Active-Link Penalty

Select the feasible set with the greatest total penalized score, not just a
greedy sequence of individually attractive connections.

For candidate connection `e`, let `score_e` be its benefit, `z_e` its binary
selection variable, and `lambda >= 0` the active-link penalty:

```text
maximize sum((score_e - lambda) * z_e)

subject to:
    z_e is 0 or 1
    z_e <= available_e
    for every terminal: sum(z_e over connections using that terminal) <= 1
```

The terminal constraint applies across both endpoint roles. A terminal cannot
transmit on one selected connection and receive on a different selected
connection at the same time.

The penalty is charged for every active connection in the selected set, not
only newly activated connections. This expresses the requested behavior:
"If a link is good but not that good, don't bother activating it."

A connection with a score below the penalty is not worth selecting under this
objective. Terminals may remain unused, and selecting no connections is allowed.
Tie handling at zero net benefit remains an implementation choice.

Use an exact optimization approach appropriate to the graph and constraints;
the particular algorithm or library has not yet been selected. The penalty must
use units compatible with the score. Its value and configuration location are
still open, since the requested policy model initially stores only the policy.

## 7. Continuous Force Evaluation Per Satellite

The explicitly preferred implementation is to loop over satellites and compute
each satellite's force separately. Returning both endpoint forces from one
pairwise function is not required.

Conceptually, at each right-hand-side (RHS) evaluation:

```text
Read the full current solver-stage state.
For each satellite:
    Initialize its total laser force to zero.
    For each active connection involving this satellite:
        Compute this connection's force on this satellite.
        Add the contribution only to this satellite's force.
    Assemble its derivatives with gravity and other existing effects.
Return the derivatives of the combined spacecraft state.
```

All satellite evaluations read the same stage snapshot. Do not advance one
satellite's position or velocity before evaluating another satellite. If the
loop is parallelized, each worker writes only its own derivative/force output
and reads a consistent shared state and active selection.

For the intended two-ended interaction, evaluate each endpoint's force with
the other endpoint as its partner. This is not double counting because each
contribution enters only that endpoint's equations. Do not also add the partner's
force inside each call, which would count it twice.

The original proposed interface was `laserForceEffector(TX, RX, linkModel)`:
calculate the acceleration on `RX`, then convert it to force using `RX` mass.
An endpoint-oriented interface, such as `force_on_endpoint(current, other, ...)`,
was suggested to avoid implying that the physical transmitter changes because
of loop order. The computational choice is settled; endpoint naming and the
precise two-ended optical force law are not.

Equal and opposite forces follow only if prescribed by that physical model;
they are not assumed merely because both endpoints are evaluated. Each
satellite uses its own mass: equal forces need not produce equal accelerations.
If the physical formula first yields force, use `acceleration = force / mass`;
if it yields acceleration, use `force = mass * acceleration` consistently.

Treat both endpoint evaluations as one selected physical connection when they
describe one interaction. Count its terminal occupancy and active-link penalty
once. The initial TX-to-RX key convention must be reconciled with this two-ended
interpretation before implementation; do not introduce a duplicate reverse
connection just to calculate the second endpoint's force.

## 8. Scheduling Every Accepted Step

Run scheduling once before the first integration step and then after every
accepted solver step. During a step, keep the selected connections fixed but
recompute their forces at every internal stage using the current stage geometry.

The intended sequence is:

1. Construct spacecraft, link parameters, registered candidates, and policy.
2. Initialize every candidate as unavailable, inactive, and zero-scored.
3. Evaluate availability at the initial state.
4. Score available candidates and solve the penalized selection problem.
5. Update active flags to match the selected set; unselected links are inactive.
6. Integrate the coupled state using continuous per-satellite laser forces.
7. After an accepted step, repeat availability, scoring, and selection.
8. Continue until simulation termination.

Internal Runge-Kutta trial stages are not accepted steps. Do not change the
schedule or commit histories during RHS evaluations; rejected steps must not
leave scheduling changes behind. When a callback changes the active selection,
use the solver's appropriate derivative/cache refresh semantics even if the
callback does not directly change positions or velocities.

Unlike the old `impulse_cb`, the new design does not apply an endpoint velocity
kick for the elapsed interval. Laser acceleration participates in the same
integration as gravity, affecting both velocity and position during that interval.

### Initial Timing Choice

The encounters discussed last approximately one hour. Start with the existing
10 s maximum orbit step:

```julia
dt_max_orbit = 10.0
```

This is a maximum, not a fixed duration. The adaptive solver can take smaller
steps. A one-hour interval contains approximately 360 steps if they are all
10 s long, with more scheduling opportunities when steps are smaller.

Ten seconds is a starting choice, not a demonstrated accuracy guarantee. Compare
against 5 s at the same solver tolerances, checking orbital effects and switching
times. Inspect actual accepted-step intervals rather than saved output times;
if the solver already takes steps below both caps, the comparison does not test
different scheduling resolutions.

Endpoint-only checks may act on a persistent boundary crossing up to roughly
one step late. If exact range/visibility boundaries are required, use event
handling. A fixed periodic scheduling cadence was discussed as an alternative,
but it is not the selected starting approach.

## 9. Details Still To Decide

- Defaults and validation rules for `n_terminal` and `health_index`, including
  the health threshold and its configuration location.
- Whether graph satellite identifiers are persistent IDs or vector indices,
  and how they map to the propagated state.
- Candidate registration rules and the representation of forbidden pairs.
- The final endpoint-direction convention and physical force law, including
  any endpoint-specific parameters.
- Exact score formulas, constellation-wide aggregation, and penalty value.
- The exact optimizer and deterministic tie handling.
- How laser-effect diagnostics are integrated and saved without restoring
  velocity kicks or accumulating rejected trial-stage contributions.

These are unresolved implementation details, not permission to add unrelated
features or to replace the agreed computational flow.

## 10. Focused Verification Goals

When implementing, verify that:

- Registered candidates start unavailable and inactive with zero scores.
- Range, health, and forbidden-pair checks prevent invalid activation.
- No terminal participates in more than one selected connection.
- A satellite can use multiple distinct terminals when permitted.
- Small graph examples reach the true maximum penalized total, including a
  case where greedy selection fails and a case where no link is worthwhile.
- Each endpoint receives exactly its own force contribution from the shared
  stage state, with correct direction and mass dependence.
- Reordering per-satellite evaluation does not change the physical result,
  apart from normal floating-point differences.
- Rejected steps do not commit schedule changes or diagnostic accumulation.
- Active-selection changes refresh the force seen by subsequent integration.
- The 10 s versus 5 s comparison meets explicitly chosen outcome tolerances.

The resulting design separates discrete connection selection from continuous
coupled dynamics while retaining the requested per-satellite force calculation.