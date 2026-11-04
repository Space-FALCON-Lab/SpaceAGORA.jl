# SMA Scheduler Scoring: Prototype, v1, and 9

## 1. Main distinction

All three SMA schedulers use the same underlying orbital relationship:
**a force component along the target's velocity increases its osculating
semimajor axis (SMA).**

The difference is what they use as the ranking score:

| Implementation | Candidate score | Meaning |
|---|---|---|
| Prototype | GVE evaluated with a unit force direction | SMA-rate sensitivity per unit acceleration |
| v1 | GVE evaluated with a unit force direction | SMA-rate sensitivity per unit acceleration |
| 9 | SMA rate evaluated with the actual candidate force and target mass | Instantaneous laser-induced SMA rate |

Thus, 9 does **not** use different orbital physics. It includes the force
magnitude and mass in its score, whereas the prototype and v1 omit that
scaling when ranking helpers. Their dynamics still apply the actual force
after selection.

The equations below describe a bound, nondegenerate osculating orbit about
a central body. They give the instantaneous contribution of the candidate
laser acceleration, not a prediction of total SMA gain over a future contact.

## 2. Start with the Gauss variational equation

For a perturbing acceleration, the SMA Gauss variational equation (GVE) is:

$$
\dot{a}_{\mathrm{laser}} =
\frac{2}{n\sqrt{1-e^2}}
\left(e\sin\nu\,A_R+\frac{p}{r}A_T\right)
$$

Here:

| Symbol | Meaning |
|---|---|
| $a$ | Target's osculating semimajor axis |
| $e$ | Target's eccentricity |
| $\nu$ | Target's true anomaly |
| $r=\lVert\mathbf{r}\rVert$ | Target's distance from the central body |
| $\mu$ | Central body's gravitational parameter |
| $n=\sqrt{\mu/a^3}$ | Mean motion |
| $p=a(1-e^2)$ | Semilatus rectum |
| $A_R$, $A_T$ | Radial and transverse components of the actual laser acceleration |

The transverse direction is in the orbital plane and perpendicular to the
radial direction, oriented along orbital motion. It is not generally parallel
to the full velocity in an eccentric orbit. The normal acceleration has no
direct contribution to instantaneous SMA rate.

Capital $A_R$ and $A_T$ denote **accelerations**, not components of a unit vector.

## 3. Separate direction from acceleration magnitude

For a candidate helper $j$, write the target's laser acceleration as:

$$
\mathbf{A}_j=\frac{\mathbf{F}_j}{m}
=\alpha_j\hat{\mathbf{u}}_j,
\qquad
\alpha_j=\frac{\lVert\mathbf{F}_j\rVert}{m}
$$

Here $m$ is the target mass, $\alpha_j$ is the acceleration magnitude, and
$\hat{\mathbf{u}}_j$ is the unit force direction on the target. In the
repulsive laser-link model used here:

$$
\hat{\mathbf{u}}_j =
\frac{\mathbf{r}_{\mathrm{target}}-\mathbf{r}_{\mathrm{helper},j}}
{\lVert\mathbf{r}_{\mathrm{target}}-\mathbf{r}_{\mathrm{helper},j}\rVert}
$$

Define its dimensionless radial and transverse components as $u_{R,j}$ and
$u_{T,j}$. Substituting $A_R=\alpha_j u_{R,j}$ and
$A_T=\alpha_j u_{T,j}$ into the GVE gives:

$$
\dot{a}_{\mathrm{laser},j}=\alpha_j q_j
$$

where

$$
q_j =
\frac{2}{n\sqrt{1-e^2}}
\left(e\sin\nu\,u_{R,j}+\frac{p}{r}u_{T,j}\right)
$$

**This $q_j$ is the quantity previously called `q_direction` in this note.**
It is explanatory notation, not a source-code variable name. It is the
SMA-rate sensitivity **per unit acceleration**, not per unit force.

Its units are seconds: multiplying $q_j$ by an acceleration in
$\mathrm{m/s^2}$ gives an SMA rate in $\mathrm{m/s}$. It depends on both the
target's orbital state and the candidate force direction; it is not merely
an angle or a dimensionless alignment score.

For positive acceleration magnitude:

- $q_j>0$: the link increases instantaneous SMA.
- $q_j=0$: the link has no instantaneous SMA effect.
- $q_j<0$: the link decreases instantaneous SMA.

## 4. Why 9's vector expression is the same GVE

9 calculates the candidate SMA rate using position and velocity directly:

$$
\dot{a}_{\mathrm{laser},j}
=\frac{2a^2}{\mu}\left(\mathbf{v}\cdot\mathbf{A}_j\right)
$$

To see why this equals the GVE, let
$h=\lVert\mathbf{r}\times\mathbf{v}\rVert=\sqrt{\mu p}$ be the target's
specific angular momentum. The radial and transverse velocity components are:

$$
v_R=\frac{\mu}{h}e\sin\nu,
\qquad
v_T=\frac{h}{r}=\frac{\mu}{h}\frac{p}{r}
$$

Therefore:

$$
\begin{aligned}
\frac{2a^2}{\mu}\left(\mathbf{v}\cdot\mathbf{A}_j\right)
&=\frac{2a^2}{\mu}\left(v_R A_R+v_T A_T\right)\\
&=\frac{2a^2}{h}\left(e\sin\nu\,A_R+\frac{p}{r}A_T\right)\\
&=\frac{2}{n\sqrt{1-e^2}}
\left(e\sin\nu\,A_R+\frac{p}{r}A_T\right).
\end{aligned}
$$

The last step uses $h=na^2\sqrt{1-e^2}$. Consequently, the directional
sensitivity also has two equivalent forms:

$$
q_j =
\frac{2a^2}{\mu}\left(\mathbf{v}\cdot\hat{\mathbf{u}}_j\right)
=
\frac{2}{n\sqrt{1-e^2}}
\left(e\sin\nu\,u_{R,j}+\frac{p}{r}u_{T,j}\right)
$$

The vector form avoids explicitly computing eccentricity and true anomaly.
For a circular orbit, the GVE reduces to $q_j=2u_{T,j}/n$; the vector form
remains usable even though true anomaly relative to periapsis is undefined.

### Why does 9 use the vector form instead of the familiar GVE?

**9 does not abandon the GVE physics.** As the derivation above shows, its
vector expression evaluates the same instantaneous SMA rate.

The vector form has practical advantages:

- It uses quantities already available to the simulation: position,
  velocity, target mass, and candidate force.
- It avoids constructing the radial/transverse frame and explicitly
  calculating eccentricity and true anomaly solely for SMA scoring.
- It works naturally for circular orbits without referring to an undefined
  periapsis. The underlying GVE remains valid in its circular limit.

These are properties of the implementation, not a claim about the original
author's intent.

9 could instead project the actual acceleration $\mathbf{F}_j/m$ into radial
and transverse components and substitute them into the familiar GVE. For the
same state and gravitational parameter, it would obtain the same rate apart
from numerical differences. Changing only the equation's representation
would therefore not change the intended scoring behavior.

The important distinction is **not "GVE versus no GVE"**. It is whether the
score includes the actual acceleration magnitude. Using a unit direction
alone would produce $q_j$, not the physical rate $\dot{a}_{\mathrm{laser},j}$.

## 5. What the implementations actually score

### Prototype and v1

The prototype's `_gve_score_laser` and v1's `_gve_score`, in their SMA branches,
project a **unit force direction** into the target's orbital frame and insert
those components into the GVE. They therefore return the directional
sensitivity $q_j$, rather than the physical rate $\alpha_j q_j$.

Although their local variables are named `aR` and `aT`, those variables
contain dimensionless direction components in these scoring functions.
They are not the actual laser accelerations.

### 9

For a single selected target, `score_gve_sma!` obtains the actual candidate
force and computes:

$$
s_j=\frac{2a^2}{\mu}
\frac{\mathbf{v}\cdot\mathbf{F}_j}{m}
=\alpha_j q_j
$$

With 9's current laser-force model and nonnegative link parameters:

$$
\lVert\mathbf{F}_j\rVert=\frac{B_jP_j}{c},
\qquad
s_j=q_j\frac{B_jP_j}{c\,m}
$$

Here $B_j$ is the edge's magnification, $P_j$ is its input power, and $c$ is
the speed of light. The score has units of $\mathrm{m/s}$.

## 6. When the helper rankings agree or differ

For the **same target state, gravitational parameter, force directions, and
eligible helper set**, if every candidate has the same positive acceleration
magnitude $\alpha$, then:

$$
s_j=\alpha q_j
\quad\Longrightarrow\quad
\operatorname*{arg\,max}_j s_j
=\operatorname*{arg\,max}_j q_j
$$

Thus, the mathematical ranking agrees. Equal-score ties may still be handled
differently. This statement concerns scoring at one instant, not identical
simulation trajectories or activation times: availability rules, update
timing, constants, and numerical treatment can differ. In particular, v1's
GVE routine uses a hard-coded Earth gravitational parameter.

If candidate acceleration magnitudes differ, the rankings can differ:

| Helper | Directional sensitivity $q_j$ (s) | Acceleration magnitude $\alpha_j$ ($\mathrm{m/s^2}$) | Physical score $s_j$ ($\mathrm{m/s}$) |
|---|---:|---:|---:|
| A | 1000 | 0.00001 | 0.010 |
| B | 800 | 0.00002 | 0.016 |

These are illustrative values, not the five-orbit simulation outputs.
Direction-only scoring prefers A. Physical-rate scoring prefers B because
its larger acceleration outweighs its less favorable direction.

This example describes the effect of unequal link strengths on the scoring
rules; it does not imply that each older implementation already exposes
per-helper power settings. v1 stores its force parameters on the link model.

For one fixed target, its mass is common to all candidates, so mass alone
cannot change their ranking. Different target masses can affect relative
scores when comparing edges involving different targets. Increasing a
positive force magnitude cannot turn a negative $q_j$ into a positive score.

## 7. Scoring is separate from matching in 9

The discussion above compares individual links for one target. 9 also has
a graph-level selection step:

- With `target_idx=i`, an edge is scored only for target `i`. Edges not
  incident to that target receive zero score.
- With `target_idx=nothing`, the policy sums the SMA-rate contributions
  of both endpoints.
- With `target_idxs` supplied, it takes precedence over `target_idx`;
  only endpoints in that selected target set contribute.
- The matcher selects a maximum-total-score set of available, positive-score
  edges, with **at most one active link per terminal**.

In the single-target example, where every candidate shares the same target
terminal, matching reduces to choosing at most one helper. This is not a
general one-link-per-spacecraft constraint: different terminals can support
different links. Only explicitly registered pairs are considered; adding
helpers does not automatically add helper-to-helper edges.

**Bottom line:** the prototype and v1 rank the GVE's directional sensitivity;
9 ranks the same physics after force/mass scaling, then resolves terminal
conflicts through matching.

## Implementation references

- Prototype direction-only score: [`2_Laser_Forces_ver2.jl`](../1_Kuang's%20Prototype%20Code/functions/2_Laser_Forces_ver2.jl)
- v1 direction-only score and force model: [`laser_link_effectors.jl`](../4_SpaceAGORA.jl-main_v1/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl)
- 9 SMA scoring: [`3_scheduling_policies.jl`](../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/3_scheduling_policies.jl)
- 9 edge force: [`laser_force_effectors.jl`](../9_SpaceAGORA.jl-main_full%20integration/src/models/dynamics/coupled/force_torque_models/laser_force_effectors.jl)
- 9 matching and scheduling: [`2_dynamic_programming_models.jl`](../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/2_dynamic_programming_models.jl)
- Five-orbit, equal-link-parameter results: [`SMA_SCHEDULER_COMPARISON.md`](SMA_SCHEDULER_COMPARISON.md)
