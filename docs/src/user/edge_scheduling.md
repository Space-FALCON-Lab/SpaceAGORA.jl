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
4. Project the window's end state on the coordinates that matter, usually the
   motion relative to the barycentre, and form the Gramian.
5. Choose the active links `σ_k ⊆ E_k`, with `|σ_k| ≤ n_max`, to maximize a
   metric of the Gramian (`schedule_edges`).

The pieces meet only through plain data: the graphs, and the matrices of the
linear model. A Cartesian, orbital-element or Kustaanheimo–Stiefel model plugs
in the same way.

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
beam = GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.05, aperture_radius_m=0.1)
link = LaserBeamLink(power_w=10_000.0, magnification=100.0, range_m=200e3, beam=beam,
    occultation=(; min_ray_altitude_m=100e3))
fore = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=deg2rad(80.0))
aft = OpticalTerminal(boresight=[0.0, -1.0, 0.0], field_of_regard_rad=deg2rad(80.0))
feasibility = EdgeFeasibilityModel(link, [[fore, aft] for _ in spacecraft];
    max_range=MinimumForceRange(min_force_n=1e-3),
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
- `nothing`, the default: only the link's `range_m`.

`power_limited_range_m(rule, link, power_w)` evaluates the rule. Each pair's
power is the link's `power_w` unless `power_w` gives one per pair, as a vector
or a function `(i, j) -> P`.

### Terminals and the field of regard

Each spacecraft has a list of `OpticalTerminal`s. A terminal steers within
`field_of_regard_rad` of its `boresight`, given in one of three frames:

- `:rtn`, the default: the spacecraft's radial, along-track and cross-track
  axes, for a spacecraft that holds a local-orbital attitude;
- `:body`: through the sample's attitude;
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

For stacked Cartesian states `[r_i; v_i]`, two helpers build the model:

- `edge_input_matrix(graph, masses)` gives the continuous input matrix of the
  link forces at the graph's geometry. Each link pushes its ends apart along
  the line of sight, as the laser links do in a run.
- `zoh_discretize(Ac, Bc, dt)` holds the inputs over a step.

```julia
steps = map(zip(midpoints, graphs)) do (t, graph)
    Ac = cartesian_jacobian(solution(t))          # the model's own Jacobian
    A, B = zoh_discretize(Ac, edge_input_matrix(graph, masses), dt)
    EdgeLinearStep(A, B, graph.pairs)
end
transform = barycentric_transform(masses)
model = EdgeLinearModel(steps, transform.T[transform.relative_rows, :])
```

A model in other coordinates supplies its own `A_k` and per-edge `B_k`. Its
projection carries its output map to Cartesian coordinates, followed by the
relative rows: `projection = T[relative_rows, :] * C`. A physical output
`δy = C δz_rel`, such as edge ranges or an RTN formation state, gives
`projection = C * T[relative_rows, :]`, and the Gramian is then the output
Gramian.

## What internal forces cannot reach

Each link force acts equally and oppositely on its two ends, so the sum of the
forces is zero and the barycentre cannot be accelerated. `barycentric_transform`
separates the barycentre rows from the relative ones. The barycentre rows of
every `B_k` are zero, and a projection onto the relative rows drops a mode that
no schedule can move.

The forces also act along the line of centres, so in free space they conserve
the constellation's angular momentum about its barycentre, and its rigid
rotation is not controllable either. In orbit, the gravity of the central body
differs between the spacecraft, which couples these modes weakly to the shape:
over a step held at constant force, the barycentre rows of a discretized `B_k`
are small but not zero. Projecting on the relative rows keeps the analysis on
the shape, whatever that coupling.

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
- `DirectionEnergyMetric(d)`: `dᵀ W d`, the reach along one chosen deformation.

`gramian_metrics(W; epsilon, direction)` returns all four at once, with the
eigenvalues, the eigenvectors and the dominant mode. `ε` is in the units of
`W`, so take it as a small share of the Gramian's mean eigenvalue.

## Schedules

`schedule_edges(model, strategy; max_active, graphs)` returns an
`EdgeSchedule`. It holds the active links of each step, the Gramian of the
schedule, and the terminals that carry each link. The strategies are:

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
many. Without `graphs`, a Boolean `feasible` matrix (edges by steps) gives the
feasible edges. `duty_cycle` and `previous` apply a duty cycle across window
boundaries.

## Minimum input energy

`minimum_input_energy(model, schedule, target; initial_state)` returns the
least `Σ |δu|²` that moves the projected end state to `target`. When an
`initial_state` is given, the inputs also correct its free drift. The result
holds the energy, the inputs by step and link, and the residual. A target
outside the reachable subspace is reported as not reachable, with its
least-squares inputs.

## Limits

- The Gramian treats each input as a signed perturbation about a nominal
  power. A laser that only pushes, from zero power, is a positive-input system,
  for which a full-rank Gramian is necessary for reachability but not
  sufficient.
- The inputs are in the units of the model's columns: one newton held over a
  step, for `edge_input_matrix` with `zoh_discretize`.
- The greedy schedules carry no optimality guarantee under terminal and
  duty-cycle limits. `log det(W + εI)` is monotone and submodular in the
  activations, and that is the property greedy bounds rest on.
- The graph's flags are evaluated at the sample times only. A link that closes
  or opens between samples is not seen.

`examples/constellation_edge_scheduling.jl` runs the pieces end to end on four
spacecraft. It compares a constant set and the three greedy schedules by their
Gramians and by the least energy of one reconfiguration.
