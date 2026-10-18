# Inclination Scheduling: Prototype, v1, and 9

## 1. What is actually implemented?

**All three now support inclination scheduling.** The original 9 snapshot
lacked it; the current workspace adds `:gve_inc` using its existing
physical-rate scoring and terminal matching pipeline.

| Code | Inclination policy | Formula used for ranking | Includes actual force magnitude and target mass? |
|---|---|---|---|
| Kuang prototype | `:gve_inc` | Inclination GVE with a unit force direction | No |
| 4 / v1 | `:gve_inc` | Inclination GVE with a unit force direction | No |
| 9 integration | `:gve_inc` (newly added) | Vector inclination rate, equivalent to the GVE away from equatorial boundaries | Yes |

The prototype and v1 use the same mathematical directional scoring rule for
inclination. This does not imply identical trajectories or event timings,
because their integration paths, scheduler update timing, and numerical
details can differ.

9 accepts `:gve_inc` alongside `:gve_sma`, `:gve_eccentricity`, and
`:gve_sma_first_index`. Like its SMA policy, it supports `target_idx`,
`target_idx=nothing` (both endpoints), and an overriding `target_idxs` set.
It scores each selected endpoint with the actual candidate force and mass,
sums contributions per edge, then selects a maximum-total-score terminal
matching of available positive-score edges.

This note documents the inspected source and the relevant mathematics.
The three-way inclination experiment is recorded separately in the
[50-helper, 50-orbit inclination results](INCLINATION_COMPARISON_50_HELPERS_50_ORBITS.md).

## 2. The physical inclination GVE

For a nondegenerate, nonequatorial osculating orbit, the instantaneous
inclination change caused by a perturbing acceleration is:

$$
\boxed{
\dot{i}=\frac{r\cos u}{h}A_N
=\frac{r\cos u}{na^2\sqrt{1-e^2}}A_N
}
$$

| Symbol | Meaning |
|---|---|
| $i$ | Inclination relative to the reference plane, between 0 and $\pi$ |
| $r=\lVert\mathbf{r}\rVert$ | Distance from the central body |
| $u=\omega+\nu$ | Argument of latitude, measured from the ascending node |
| $\mathbf{h}=\mathbf{r}\times\mathbf{v}$ | Specific angular momentum vector |
| $h=\lVert\mathbf{h}\rVert$ | Specific angular momentum magnitude |
| $a$, $e$ | Osculating semimajor axis and eccentricity |
| $n=\sqrt{\mu/a^3}$ | Mean motion for a bound orbit |
| $\hat{\mathbf{n}}=\mathbf{h}/h$ | Unit normal to the target's orbital plane |
| $A_N=\mathbf{A}\cdot\hat{\mathbf{n}}$ | Normal component of the actual perturbing acceleration |

The relation $h=na^2\sqrt{1-e^2}$ explains the two equivalent coefficients.
The expression involving $h$ does not require explicitly computing $a$,
$e$, or $n$.

Only the acceleration **normal to the orbital plane** contributes directly
to instantaneous inclination rate. Radial and transverse accelerations do
not directly rotate the orbital plane.

The factor $\cos u$ matters:

- At the ascending node ($u=0$), positive $A_N$ increases inclination.
- At the descending node ($u=\pi$), negative $A_N$ increases inclination.
- At $u=\pi/2$ or $3\pi/2$, the instantaneous inclination contribution is zero.
  A normal force there can instead change the node direction.

Thus, "choose the largest positive normal force" is not generally the same
as "choose the largest positive inclination rate."

## 3. What the prototype and v1 actually score

For candidate helper $j$, both implementations construct the unit force
direction on the target:

$$
\hat{\mathbf{f}}_j=
\frac{\mathbf{r}_{\mathrm{target}}-\mathbf{r}_{\mathrm{helper},j}}
{\lVert\mathbf{r}_{\mathrm{target}}-\mathbf{r}_{\mathrm{helper},j}\rVert}
$$

They project it onto the target's orbital normal:

$$
f_{N,j}=\hat{\mathbf{f}}_j\cdot\hat{\mathbf{n}}.
$$

Their local variable named `aN` holds this **dimensionless direction
component**, not an actual acceleration. The returned inclination score is:

$$
\boxed{
q_{i,j}=\frac{r\cos u}{h}f_{N,j}
}
$$

$q_{i,j}$ is notation used in this note, not a source-code variable name.
It is the inclination-rate sensitivity **per unit acceleration magnitude**.
With radians treated as dimensionless, it has units of $\mathrm{s/m}$.

To obtain the physical rate, write
$\mathbf{A}_j=\alpha_j\hat{\mathbf{f}}_j$, where
$\alpha_j=\lVert\mathbf{F}_j\rVert/m$:

$$
\boxed{
\dot{i}_j=\alpha_j q_{i,j}
=\frac{\lVert\mathbf{F}_j\rVert}{m}q_{i,j}
}
$$

The result is in radians per second. It is not in degrees per second.
For a sensitivity per unit **force**, the coefficient would instead be
$q_{i,j}/m$.

Both older schedulers choose the eligible helper with the largest
**strictly positive** score, or choose no helper if none is positive.
Their force models subsequently apply the actual force in the dynamics.
They maximize signed inclination increase, not absolute plane rotation,
inclination decrease, or tracking toward a prescribed inclination.

### Selection and implementation details

| Detail | Prototype | v1 |
|---|---|---|
| Scorer | `_gve_score_laser`, `:gve_inc` branch | `_gve_score`, `:gve_inc` branch |
| Selection | `_gve_select_best_cavity` picks one target-helper cavity | `update_laser_link_schedule!` picks one helper for its target |
| Availability | Configured range and optional line-of-sight/clearance checks | Configured range |
| Update path in the existing comparison runners | Re-evaluated during RHS force evaluation | Re-evaluated at accepted solver steps |
| Gravitational parameter used by scorer | Supplied `mu` | Hard-coded `_MU_EARTH_GVE` |

For a shared target state, gravitational parameter, and force direction,
their mathematical inclination scoring formulas agree. Their orbital-element
reconstruction includes numerical guards and eccentricity clamping, so this
is not a claim of exact floating-point equality in every state.

Equal positive acceleration magnitudes preserve the ranking between
$q_{i,j}$ and $\dot{i}_j$. Unequal candidate acceleration magnitudes can
change that ranking. This is the same scaling distinction discussed for
[SMA scoring](SMA_SCHEDULER_SCORING.md); 9 now includes that physical scaling
for inclination as well.

## 4. Equivalent vector form and 9's implementation

An inclination-rate scorer can also be expressed using Cartesian vectors.
Let $\hat{\mathbf{k}}$ be the fixed unit normal to the reference plane
(the positive inertial z-axis here). Then:

$$
\cos i=\frac{\hat{\mathbf{k}}\cdot\mathbf{h}}{h},
\qquad
\dot{\mathbf{h}}=\mathbf{r}\times\mathbf{A}.
$$

The second expression gives the candidate perturbing acceleration's
contribution to angular-momentum change; central gravity contributes no
torque about the central body.

Differentiating $\cos i$ gives:

$$
-\sin i\,\dot{i}
=
\frac{\hat{\mathbf{k}}\cdot\dot{\mathbf{h}}}{h}
-
\frac{\cos i}{h}\left(\hat{\mathbf{n}}\cdot\dot{\mathbf{h}}\right).
$$

For $h>0$ and $0<i<\pi$:

$$
\boxed{
\dot{i}=
\frac{
\cos i\,\left(\hat{\mathbf{n}}\cdot(\mathbf{r}\times\mathbf{A})\right)
-\hat{\mathbf{k}}\cdot(\mathbf{r}\times\mathbf{A})
}{h\sin i}
}
$$

This is equivalent to the physical inclination GVE, not a different
objective. To see this, let $\hat{\mathbf{t}}$ be the transverse unit vector.
Then:

$$
\mathbf{r}\times\mathbf{A}
=rA_T\hat{\mathbf{n}}-rA_N\hat{\mathbf{t}},
\qquad
\hat{\mathbf{k}}\cdot\hat{\mathbf{t}}=\sin i\cos u.
$$

Substitution cancels the $A_T$ terms and leaves
$\dot{i}=r\cos u\,A_N/h$.

Using the actual acceleration $\mathbf{F}_j/m$ in this vector expression
would give the physical rate. Using only $\hat{\mathbf{f}}_j$ would give
the directional sensitivity $q_{i,j}$ instead.

9 evaluates the same derivative through an `atan` representation rather
than dividing explicitly by $\sin i$. Define:

$$
\rho=\sqrt{h_x^2+h_y^2},\qquad i=\operatorname{atan2}(\rho,h_z).
$$

For $\rho>0$:

$$
\dot{\rho}=\frac{h_x\dot{h}_x+h_y\dot{h}_y}{\rho},
\qquad
\dot{i}=\frac{h_z\dot{\rho}-\rho\dot{h}_z}{h^2}.
$$

The implementation uses normalized factors to evaluate this expression
without forming $h^2$. It requires positive finite mass and finite nonzero
angular momentum, and reports invalid input as an error. The rate is in
radians/s. Its GVE equivalence, scaling, boundary behavior, target filtering,
and matching are covered by dedicated tests.

## 5. Circular and equatorial orbits need different care

### Circular but inclined

For $e=0$, periapsis and therefore $\omega$ and $\nu$ separately are
undefined. However, the ascending node and argument of latitude $u$ remain
defined when the orbit is not equatorial.

Both older scorers handle this case by computing $u$ from the node direction
and position rather than relying on a periapsis direction.

### Exactly equatorial

At $i=0$ or $\pi$, the ascending node is undefined. Both older scorers use
an equatorial fallback:

```julia
atan(tgt_pos[2], tgt_pos[1])
```

This substitutes true longitude from the x-axis for argument of latitude.
The `:gve_inc` branch does not reject the equatorial case, unlike the nearby
RAAN branch's small-inclination guard.

**That fallback avoids an undefined node calculation but does not remove
the mathematical ambiguity of a signed inclination rate at an exactly
equatorial orbit.** Its score should not be interpreted as a generally
valid physical inclination derivative for an arbitrary out-of-plane force
at that boundary.

In particular, from an exactly prograde equatorial orbit, either sign of a
nonzero normal acceleration can initially tilt the plane away from the
reference plane. The unsigned inclination magnitude has a one-sided growth
rate; a longitude-dependent signed score need not represent it correctly.

The vector formula above also divides by $\sin i$, so it is singular at
equatorial inclination and can be poorly conditioned nearby. Unlike the
SMA vector expression, it does **not** automatically solve the singularity.
9 explicitly handles $\rho=0$ using the forward one-sided derivative:

$$
\dot{i}_{+}=\operatorname{sign}(h_z)
\frac{\sqrt{\dot{h}_x^2+\dot{h}_y^2}}{h}.
$$

Thus either normal force sign can increase inclination from zero; at pi
either plane tilt decreases inclination. A coplanar force has zero boundary
rate. There is no arbitrary node direction or small-inclination cutoff.
Away from the exact boundary, the signed derivative is used. This boundary
behavior intentionally differs from the older longitude fallback.

### Relevance to the previous simulations

The earlier 20-helper and 50-helper simulations used coplanar equatorial
helpers and target. Their separation vectors, and hence laser forces, lie
in the common orbital plane. In the ideal geometry:

$$
A_N=0
\quad\Longrightarrow\quad
\dot{i}=0.
$$

That geometry cannot provide a meaningful inclination-raising comparison.
An inclination experiment would need out-of-plane link geometry, preferably
with a nonzero initial target inclination to avoid conflating scheduling
behavior with equatorial singularities. The existing simulations used SMA
scheduling; they are not inclination test results.

## 6. Bottom line

- **Prototype and v1:** implemented inclination GVE schedulers, both ranking
  unit-direction sensitivity rather than force/mass-scaled physical rate.
- **9:** physical inclination-rate scheduler with actual force/mass scaling,
  shared target filtering with SMA, and terminal matching.
- **GVE versus vector form:** equivalent physical rate expressions exist,
  but neither should be confused with the unscaled directional score.
- **Equatorial behavior:** needs explicit care; switching to a vector
  expression alone is not sufficient.

## 7. Inclination support across the other snapshots

The workspace has 11 numbered snapshot directories because there are three
variants prefixed with `8`. This inventory records source implementation,
not a certification that each snapshot's complete application runs.

| Snapshot | Inclination scheduling source present? |
|---|---|
| [1 / prototype](../1_Kuang's%20Prototype%20Code/) | Yes, `:gve_inc` |
| [2 / SpaceAGORA](../2_SpaceAGORA.jl/) | Yes, `:gve_inc` |
| [3 / blank](../3_SpaceAGORA.jl-main_blank/) | No inclination scheduler found |
| [4 / v1](../4_SpaceAGORA.jl-main_v1/) | Yes, `:gve_inc` |
| [5 / v2](../5_SpaceAGORA.jl-main_v2/) | Yes, `:gve_inc` |
| [6 / v3](../6_SpaceAGORA.jl-main_v3/) | Yes, `:gve_inc` |
| [7 / v4](../7_SpaceAGORA.jl-main_v4/) | Yes, `:gve_inc` |
| [8 / v5 single target](../8_SpaceAGORA.jl-main_v5_single_target/) | No inclination policy |
| [8 / v5 indexed single target](../8_SpaceAGORA.jl-main_v5_single_target_with_index_edge_for_speed/) | No inclination policy |
| [8 / v6 indexed multi-target](../8_SpaceAGORA.jl-main_v6_multi_target_with_index_edge_for_speed/) | No inclination policy |
| [9 / full integration](../9_SpaceAGORA.jl-main_full%20integration/) | Yes, `:gve_inc` added in this workspace |

### How 7 implements it

In [7's laser-link source](../7_SpaceAGORA.jl-main_v4/I_src/1.3_dynamics/coupled/force_torque_models/laser_link_effectors.jl),
`choose_active_links!` implements this sequence:

1. Read the shared schedule setting from the first registered laser model.
2. Collect registered pairs within their individual `range_m`.
3. Treat each pair `(sat_i, sat_j)` as `(receiver, emitter)` for scoring.
4. Compute the receiver's directional GVE score
   $q_i=(r\cos u/h)(\hat{\mathbf{f}}\cdot\hat{\mathbf{n}})$.
5. Discard nonpositive scores and sort the remainder from highest to lowest.
6. Greedily select links whose spacecraft have not already been assigned.

It reselects from scratch, rather than retaining a previous link solely
because it remains in range. The scheduler callback runs at initialization
and every accepted ODE step.

This is **satellite-disjoint greedy selection**, not maximum-total-score
terminal matching. Multiple independent links can be active, but a spacecraft
can appear in at most one selected pair. Greedy selection is not guaranteed
to maximize the sum of scores across a general network. If every candidate
shares one target, it reduces to choosing the highest positive-scoring link.

The score uses a unit force direction and omits power, magnification, and
mass. Actual force enters a separate `laser_impulse_callback`: it computes
force magnitude $\eta\beta M P/c$, divides by the supplied `mass_kg`,
multiplies by the elapsed timestep, and applies opposite velocity increments
to receiver and emitter. The force effector's `wrench` itself returns zero;
this implementation applies laser action as discrete velocity kicks rather
than a continuous laser acceleration in the ODE RHS.

7 uses the same longitude fallback at equatorial inclination, so the
boundary caveat in section 5 also applies.

## Source references

- [Prototype scoring, selection, and force direction](../1_Kuang's%20Prototype%20Code/functions/2_Laser_Forces_ver2.jl)
- [v1 scoring, policy selection, and force model](../4_SpaceAGORA.jl-main_v1/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl)
- [9 supported policies and dispatcher](../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/2_dynamic_programming_models.jl)
- [9 implemented orbital-element scorers](../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/3_scheduling_policies.jl)
- [SMA scoring explanation](SMA_SCHEDULER_SCORING.md)
- [50-helper SMA simulation comparison](SMA_SCHEDULER_COMPARISON_50_HELPERS_50_ORBITS.md)
