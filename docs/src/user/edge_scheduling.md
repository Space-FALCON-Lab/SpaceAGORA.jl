# Edge Scheduling

This page covers planning which laser links of a constellation to activate over
a horizon so that its formation is easy to reshape: the feasible-edge graph of
each step, a linear model of the window, the finite-horizon controllability
Gramian, and schedules that make it large.

What to read next:

- [Interlinks](interlinks.md)
- [Constellations](constellations.md)

## From a propagated constellation to a schedule

A run's interlink scheduler picks links step by step from what each link does
now (see [Interlinks](interlinks.md#Scheduling)). Edge scheduling plans over a
window instead. The pieces are:

1. Propagate the constellation, for instance with
   `run_simulation(args; return_solution=true)`.
2. At each planning step `k`, evaluate the feasible-edge graph `G_k`: which
   spacecraft pairs could carry a link (`feasible_edge_graph`).
3. Build a linear time-varying model of the window,
   `δx_{k+1} = A_k δx_k + Σ_e B_k^e δu_{k,e}`, with one input block per
   candidate link, in any state representation (`EdgeLinearModel`).
4. Project the window's end state on the coordinates that matter: the
   formation's shape, in one unit (see
   [What the links cannot reach](#What-the-links-cannot-reach)).
5. Choose the active links `σ_k ⊆ E_k`, with `|σ_k| ≤ n_max`, to maximize a
   metric of the Gramian (`schedule_edges`).

The pieces meet only through plain data: the graphs, and the matrices of the
linear model. A Cartesian, orbital-element or Kustaanheimo–Stiefel model plugs
in the same way. In a model whose spacecraft keep their own time variables,
such as Kustaanheimo–Stiefel fictitious time, each step must span one common
physical-time interval for every spacecraft, so that both ends of a link are
forced together.

## The feasible-edge graph

An `EdgeFeasibilityModel` holds the rule. An edge `(i, j)` is feasible when:

- **Line of sight.** With the link's `occultation`, the segment between the two
  spacecraft clears the planet's sphere, raised by `min_ray_altitude_m`, as for
  the link in a run.
- **Range against power.** `0 < ρ_ij ≤ min(range_m, ρ_max(P_ij))`. `P_ij` is the
  pair's transmitted power, and `ρ_max` is the `max_range` rule.
- **Field of regard at both ends.** Some terminal of `i` sees the line of sight
  to `j` within its field of regard, and some terminal of `j` sees the line of
  sight back to `i`.
- **Terminal and duty-cycle limits.** Links active together need distinct
  terminals at each end, and a `DutyCycleLimit` caps how often links may be
  active.

```julia
beam = GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.2, aperture_radius_m=0.2)
link = LaserBeamLink(power_w=10_000.0, magnification=100.0, range_m=200e3, beam=beam,
    occultation=(; min_ray_altitude_m=100e3))
fore = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=deg2rad(70.0))
aft = OpticalTerminal(boresight=[0.0, -1.0, 0.0], field_of_regard_rad=deg2rad(70.0))
feasibility = EdgeFeasibilityModel(link, [[fore, aft] for _ in spacecraft];
    max_range=MinimumForceRange(min_force_n=2e-3),
    duty_cycle=DutyCycleLimit(on_steps=4, window_steps=6))

solution = run_simulation(args; return_solution=true)
graphs = feasible_edge_graphs(feasibility, solution, times; planet=planet)
feasible_pairs(graphs[1])
```

A `FeasibleEdgeGraph` keeps each constraint's flag for every candidate pair
(`line_of_sight`, `in_range`, `field_of_regard`, `duty_cycle`), the separation
and its limit, the unit line of sight `direction`, and the terminals that can
point along it. A single sample time takes `feasible_edge_graph(feasibility,
samples; context=LinkContext(time_s=t, planet=planet))`, with one `StateSample`
per spacecraft or a run state.

The link is a `LaserBeamLink` or a `LaserMomentumLink`. Its own `occultation`
and `range_m` apply. The field of regard comes from the terminals, so the link
carries no terminal pointing or mounted recoil. `pairs` defaults to every pair
(`all_pairs`). The `ring_pairs` and `grid_pairs` lists of a constellation work
too.

### Range against power

`max_range` gives `ρ_max(P)`:

- `MinimumForceRange(min_force_n=F)`: the largest separation at which the
  link still pushes with at least `F` newtons. With a `GaussianBeam`, the force
  falls with the beam's capture fraction at the far aperture, and `ρ_max`
  follows in closed form. Without a beam, the force does not fall with range.
- Any function of the power in watts that returns metres.
- A positive number of metres: a cap that does not depend on the power.
- `nothing`, the default: only the link's `range_m`.

`power_limited_range_m(rule, link, power_w)` evaluates the rule. Each pair's
power is the link's `power_w` unless `power_w` gives one per pair, as a vector
or a function `(i, j) -> P`.

### Terminals and the field of regard

Each spacecraft has a list of `OpticalTerminal`s. A terminal steers within
`field_of_regard_rad` of its `boresight`, given in one of three frames:

- `:rtn`, the default: the spacecraft's radial, along-track and cross-track
  axes, for a spacecraft that holds a local-orbital attitude;
- `:body`: through the sample's attitude, as the link's own terminal pointing
  reads mounted hardware;
- `:inertial`.

A terminal carries one link at a time. `terminal_assignment(graph, active)`
finds distinct terminals for a set of active links, or `nothing` when none
exist. `admissible_edge_set(graph, active; max_active, previous)` checks a whole
step: line of sight, range, field of regard, the count limit, the terminals and
the duty cycle.

### Duty cycle

`DutyCycleLimit(on_steps=n, window_steps=w)` allows at most `n` active
link-steps in any `w` consecutive steps. With `scope=:spacecraft` (the default),
the count is per spacecraft, and each of its active links counts once per step.
With `scope=:edge`, the count is per link. A graph evaluated with `previous`
(the active sets of the steps before, oldest first) marks an edge infeasible
when activating it alone would exceed the limit. The scheduler applies the
limit to its whole schedule.

## The linear model of a window

An `EdgeLinearStep(A, B, edges)` is one step: the state transition `A_k` and
the input columns of every candidate edge. By default each edge has one column,
its force magnitude in newtons. Edges with more inputs give their column ranges
with `columns`, or the step takes one block per edge. An
`EdgeLinearModel(steps, projection)` is the window, and `projection` maps the
state at its end to the coordinates whose controllability is assessed.

For stacked Cartesian states `[r_i; v_i]`, these helpers build the model:

- `edge_input_matrix(graph, masses)` gives the continuous input matrix of the
  link forces at the graph's geometry. A positive input pushes the ends apart
  along the line of sight, as the laser links do in a run.
- `link_force_jacobian(graph, masses, forces)` gives the state Jacobian of
  nominal link forces, for a model linearized about a positive force (see
  [Signed and push-only inputs](#Signed-and-push-only-inputs)).
- `zoh_discretize(Ac, Bc, dt)` holds the inputs over a step.

The dynamics' own Jacobian is the model's to supply. The example
`examples/constellation_edge_scheduling.jl` builds it for point-mass gravity in
`cartesian_step`:

```julia
steps = map(zip(midpoints, graphs)) do (t, graph)
    Ac = gravity_jacobian(solution(t))           # your model's Jacobian
    A, B = zoh_discretize(Ac, edge_input_matrix(graph, masses), dt)
    EdgeLinearStep(A, B, graph.pairs)
end
model = EdgeLinearModel(steps, range_output_matrix(solution(t_end), pairs; rate_scale=1 / n))
```

A model in other coordinates supplies its own `A_k` and per-edge `B_k`. Its
projection carries its output map to Cartesian coordinates, followed by the
Cartesian output: `projection = C_out * C_cart`. A physical output of the
relative motion, `δy = C δz_rel`, gives `projection = C * T[relative_rows, :]`
with `T` from `barycentric_transform`, and the Gramian is then the output
Gramian.

## What the links cannot reach

Each link force acts equally and oppositely on its two ends along the line of
centers. Two conserved quantities follow:

- **The barycenter.** The forces sum to zero, so the continuous input matrix
  `B_c` has no barycenter component. `barycentric_transform` separates the
  barycenter rows from the relative ones.
- **The rigid rotation.** The forces also leave the total angular momentum
  about the planet's center unchanged, and in free space the angular momentum
  about the barycenter too. So the formation's rigid rotation cannot be driven
  in free space.

In orbit, unequal gravity couples both weakly to the shape. Over a step held at
constant force, the barycenter rows of a discretized `B_k` are small but not
zero, and the rotation about the barycenter can change only by trading angular
momentum with the barycenter's orbit, which is slow over a window short
against the orbital period.

The relative rows `T[relative_rows, :]` therefore drop the barycenter but keep
the rotation. On them, `λ_min` measures the weakly reachable rotation rather
than the shape, and so does `log det` with a small `ε`. For those two metrics,
project on a rotation-free output: `range_output_matrix(states, pairs)` gives
the pairwise ranges and range rates, which neither a translation nor a rotation
of the whole formation changes. `λ₁` and the direction energy can keep the
relative rows, since there the user chooses the direction.

Rows in different units spread the Gramian's eigenvalues by the units alone. A
projection with ranges in metres and range rates in metres per second does so
too, so scale the rates to metres: `rate_scale = 1 / n` for the mean motion
`n`, or divide velocity rows by `n` on a relative projection. Then choose `ε`
against the scaled spectrum.

## The Gramian and its metrics

`reachability_matrix(model; active)` returns `R_N = P [Φ(k₀+N, k₀+1) B_{k₀} …
B_{k₀+N−1}]`, restricted to the active links, together with a name for each
column. `controllability_gramian(model; active)` returns `W_N = R_N R_Nᵀ`. With
a physical output projection, the column of an activation is that activation's
terminal gain: the change in the output at the window's end per unit input.

The metrics are the values of `controllability_metric`:

- `LargestEigenvalueMetric()`: `λ₁`, the easiest direction. Its eigenvector is
  the dominant deformation mode.
- `SmallestEigenvalueMetric()`: `λ_min`, the hardest direction.
- `LogDetMetric(epsilon=ε)`: `log det(W + εI)`, the volume of the reachable set.
  Directions with eigenvalues below `ε` count little.
- `DirectionEnergyMetric(d)`: `dᵀ W d`, the reach along one chosen deformation.

`gramian_metrics(W; epsilon, direction)` returns all four at once, with the
eigenvalues, the eigenvectors and the dominant mode.

## Schedules

`schedule_edges(model, strategy; max_active, graphs)` returns an
`EdgeSchedule`. It holds the active links of each step, the Gramian of the
schedule, the terminals that carry each link, and `tied_rounds`. The strategies
are:

- `ConstantEdgeSet(edges=...)`: one link set for the whole window, each link
  active at every step where it is admissible. With `metric` in place of
  `edges`, the set is built from the candidates, one link per round, by that
  metric.
- `GreedyEdgeSchedule(LargestEigenvalueMetric())`: greedy `λ₁`.
- `GreedyEdgeSchedule(LogDetMetric(epsilon=ε))`: greedy volume.
- `GreedyEdgeSchedule(SmallestEigenvalueMetric(); tie_break=LogDetMetric(epsilon=ε))`:
  greedy `λ_min`. `λ_min` stays zero until every direction can be reached, so
  until then the tie break decides.

The greedy schedules add one `(step, link)` activation at a time over the
whole window, and fill every admissible slot. The built-in metrics never fall
as links are added, so the metric chooses which links are active, not how
many. When no limit binds, every strategy returns the same schedule, with every
feasible link active; strategies differ only under a binding `max_active`,
terminal or duty-cycle limit.

Two values count as equal within the roundoff of a symmetric eigensolver,
about `eps · n · λ_max`, carried into each metric. `tied_rounds` counts the
rounds in which every candidate tied on the metric, so that the tie break or
the candidate order chose. A greedy `λ_min` schedule whose `tied_rounds` equals
its number of activations has been chosen by its tie break alone, which is the
sign of a projection that keeps a direction no link can reach.

Without `graphs`, a Boolean `feasible` matrix (edges by steps) gives the
feasible edges. `duty_cycle` and `previous` apply a duty cycle across window
boundaries.

## Minimum input energy

`minimum_input_energy(model, schedule, target; initial_state)` returns the
least `Σ |δu|²` that moves the projected end state to `target`. When an
`initial_state` is given, the inputs also correct its free drift. The result
holds:

- the energy, the inputs by step and link, and the residual;
- whether the target is reachable: a target outside the reachable subspace gets
  its least-squares inputs and is reported as not reachable;
- `most_negative_input` and `nonnegative`, the sign check described next.

### Signed and push-only inputs

The inputs are signed in the sign of the model's columns. With
`edge_input_matrix`, a positive input pushes the ends apart and a negative one
pulls them together, which a laser cannot do. Which signed answer is
meaningful depends on the nominal force the model is linearized about:

- **About zero link force**, as in the example, the inputs are the link forces
  themselves. A push-only laser then needs `u ≥ 0`: the target is reachable
  when `R u = y` has a nonnegative solution, and a full-rank Gramian is
  necessary for that but not sufficient. When `nonnegative` is false, the
  energy is that of a two-way actuator, a lower bound on the push-only energy.
- **About a positive nominal force `F*`**, signed perturbations are legitimate
  within `−F* ≤ δu ≤ F_max − F*`. The model's `A_k` must then include the
  nominal force's own Jacobian, `(F*/ρ)(I − ℓ̂ℓ̂ᵀ)` on the relative position,
  which `link_force_jacobian` adds to `Ac`.

The Gramian and the energy do not depend on the sign convention of the columns;
the returned inputs do.

## Limits

- The greedy schedules carry no optimality guarantee under terminal and
  duty-cycle limits. `log det(W + εI)` is monotone and submodular in the
  activations, and that is the property greedy bounds rest on.
- Each greedy round evaluates every admissible `(step, link)` slot, with an
  eigendecomposition for the eigenvalue metrics. That suits formations of a
  few spacecraft over tens of steps, not a proliferated constellation.
- The graph's flags are evaluated at the sample times only. A link that closes
  or opens between samples is not seen.
- `link_force_jacobian` holds the force's magnitude fixed. A beam force that
  falls with range adds a radial term that it leaves out.

`examples/constellation_edge_scheduling.jl` runs the pieces end to end on four
spacecraft. It compares a constant set and the three greedy schedules by their
Gramians on the formation's ranges and range rates, and by the least two-way
energy of one reconfiguration, with its sign check.
