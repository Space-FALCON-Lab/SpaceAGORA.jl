# Interlinks

Use this page to connect spacecraft through interlinks, with a policy that
selects which links are active.

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [Parallel Execution](parallel_execution.md)

## Interlinks and constellations

An **interlink** is a connection between two spacecraft terminals that does
something while the geometry allows it: push the spacecraft apart, carry power,
or carry data. A **constellation** is a set of spacecraft and
their orbits, laid out and run together. They are different things: a
constellation needs no interlinks, and two spacecraft can be linked without
forming one.

## Links and link kinds

Each candidate link has its own kind, a value of a subtype of
`AbstractInterlink`. The kind decides when the link can close and what it does
to its two endpoints:

- `link_range_m(link)`: the maximum separation in metres. Every kind
  implements it.
- `link_available(link, a, b, ctx)`: whether the link can close. By default,
  the endpoints are apart and within `link_range_m`.
- `link_force(link, a, b, ctx)`: the inertial forces (N) on `a` and `b` while
  the link is active. By default, none.

`a` and `b` are the endpoints' state samples (`StateSample`: position,
velocity, mass, attitude and `spacecraft_index`), in the order of the link's
key: the lower `(spacecraft, terminal)` pair first, or, for a directed link, the
transmitter first. `ctx` is a `LinkContext`
holding the time (s into the run) and the run's planet.

The built-in kinds are:

- `LaserMomentumLink(; power_w, magnification, range_m)`: an open-cavity laser
  link with the idealized central force `magnification * power_w / c`, directed
  away from the other endpoint, equal and opposite on the two spacecraft.
  `LaserMomentumLink(terminal)` takes the numbers of a `LaserTerminal`.
- `LaserBeamLink(; power_w, magnification, range_m, ...)`: the same link, with
  opt-in occultation, beam divergence, pointing and mounted recoil. See
  [Laser physics](#Laser-physics).
- `PowerBeamLink(; power_w, wavelength_m, waist_radius_m, range_m)`: a power
  beam from a `LaserTerminal` to a `PowerReceiver`. See [Power links](#Power-links).
- `DataLink(; transmit_power_w, budget, range_m)`: a one-way data link, from a
  `RadioTerminal` or a virtual terminal to another, at the rate its link budget
  sets. See [Data links](#Data-links).

## Setting up interlinks

An `InterlinkModel` holds the candidate links, spacecraft pairs that may never
link, and the history of the links the scheduler selected. Endpoints are
`(spacecraft vector index, terminal number)` pairs: the spacecraft's position
in `dynamics_model.spacecraft`, not its `id`. Each terminal can take part in at
most one selected link.

```julia
model = InterlinkModel()
for partner in 2:3
    register_link!(model, spacecraft, (1, 1), (partner, 1),
        LaserMomentumLink(power_w=10_000.0, magnification=100.0, range_m=200e3))
end
args = SimulationConfiguration(...; interlink_model=model,
    scheduling_policy_model=SemimajorAxisRatePolicy())
```

Start with `examples/gve_sma_interlinks.jl`. It builds a three-spacecraft case,
registers two candidate links, and selects the link that increases the target's
instantaneous semimajor-axis rate. `examples/targeted_40sat_sma_interlinks.jl`
uses the same interface with thirty helpers and ten targets. These are analytical
atmosphere examples; they do not require native GRAM or private telemetry.

```julia
include("examples/gve_sma_interlinks.jl")
args = gve_sma_configuration(duration=120.0)
solution = run_simulation(args; return_solution=true)
history = solution.prob.p.args.interlink_model.history
```

## Terminals

A spacecraft's terminals are its `n_terminal` virtual terminals, numbered
`1:n_terminal`, followed by the terminal hardware mounted on it, in mounting
order. A `LaserTerminal` mounted with `add_actuator!` is terminal hardware; a
kind of your own opts in with `SpaceAGORA.link_terminal(::MyTerminal) = true`.
`terminal_count(spacecraft)` gives the total, and
`mounted_link_terminals(spacecraft)` lists the hardware with its numbers. A
spacecraft with only hardware sets `n_terminal=0`. Virtual terminals stay for
idealized studies, which run as before.

An endpoint's terminal is a number, or the mounted terminal itself:

```julia
bus = Link(root=true, m=100.0)
sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, n_terminal=0)
terminal = LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=50.0, magnification=100.0, range_m=200e3)
add_actuator!(sc, bus, terminal)
# ... a second spacecraft, `other`, likewise
spacecraft = [sc, other]
register_link!(model, spacecraft, (1, terminal), (2, 1), LaserMomentumLink(terminal))
```

A terminal must be able to carry its link (`terminal_supports(terminal, link)`).
A `LaserTerminal` carries `LaserMomentumLink`, and a virtual terminal carries
any kind, except a `DataLink` whose `ShannonBudget` reads the radios at both
ends. `terminal_number(spacecraft, terminal)` gives a mounted terminal's
number. Two mounted terminals with equal fields cannot be told apart, so for
those give the number.

### Directed links

A kind with `link_directed(link) = true` runs one way: from the first endpoint
given to `register_link!` (the transmitter) to the second (the receiver). Its
key keeps that order, so its hooks receive the transmitter as `a`. The opposite
direction between the same terminals is a separate link. Laser momentum links
are undirected.

### Registering many links

`register_links!(model, spacecraft, pairs, link)` registers one link per
spacecraft pair. `all_pairs(spacecraft)` gives every pair, or you can pass any
vector of `(i, j)` index pairs. For a constellation laid out by
`build_constellation` or `walker_constellation` (see
[Constellations](constellations.md)), its planes and slots give the usual
topologies:

- `same_plane_pairs(spacecraft)`: every pair within each plane.
- `ring_pairs(spacecraft)`: each slot with the next in its plane, closing the
  ring.
- `grid_pairs(spacecraft; wrap_planes=true)`: the "+grid", the rings plus each
  slot with the same slot of the next plane. Pass `wrap_planes=false` for a
  Walker star, whose first and last planes counter-rotate. `link` is a kind value, or a function
`(i, j) -> link`. Each spacecraft's links take its terminals that can carry them
in turn, starting again from the first when they run out; links that share a
terminal compete for it when the scheduler selects. Pass
`terminals=(i, j) -> (terminal_i, terminal_j)` to choose the terminals yourself.

```julia
register_links!(model, spacecraft, all_pairs(spacecraft), LaserMomentumLink(power_w=40.0))
# or, for a constellation's +grid:
register_links!(model, spacecraft, grid_pairs(spacecraft), LaserMomentumLink(power_w=40.0))
```

### Terminal hardware and power

A hardware terminal that is a link endpoint is driven by its links for the
whole run:

- **Its level.** While one of its links is active, it runs at the level the
  link sets, `link_terminal_levels(link)`. For a `LaserMomentumLink`, that is
  the link's `power_w` at both ends. Otherwise its level is 0.
- **Its power.** It draws `actuator_consumption(terminal, level)` from its
  spacecraft's power budget: `power_w / efficiency` for a `LaserTerminal`. The
  draw shows in the power columns, and in the `actuator_level` and
  `actuator_power` columns.
- **Its force.** Only the link's force applies while the link is active. The
  terminal's own recoil applies only when a controller drives a terminal that
  is no link's endpoint.
- **Controllers.** A mounted-actuation controller's levels for a link-driven
  terminal apply nothing, and setup says so once.
- **Limits.** A link that would drive its terminal past its limits (above a
  `LaserTerminal`'s `max_power_w`) is refused at setup. Build the link from the
  terminal, `LaserMomentumLink(terminal)`, to match them.

When the battery sheds tier 1 (see [Powering a Spacecraft](power_budget.md)), it
sheds laser and radio terminals. Links through shed hardware terminals then
become unavailable, and stay so until the tier clears: a shed terminal keeps
asking for its link's power, so the tier clears when the battery can supply it
again. A terminal kind's tier is `terminal_shed_tier(terminal)`: 1 for
`LaserTerminal` and `RadioTerminal`, 0 (never shed) for `PowerReceiver` and by
default. Virtual terminals draw nothing and are never shed.

## Laser physics

`LaserBeamLink` is the open-cavity laser link with four effects you can switch
on. With none of them it is `LaserMomentumLink`, bit for bit.

```julia
beam = GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.05, aperture_radius_m=0.08)
link = LaserBeamLink(power_w=10_000.0, magnification=100.0, range_m=200e3;
    occultation=(; min_ray_altitude_m=100e3),   # or true
    beam,
    pointing=TerminalPointing(field_of_regard_rad=deg2rad(30)),
    recoil=:mounted)
register_link!(model, spacecraft, (1, terminal_1), (2, terminal_2), link)
```

- **Line of sight** (`occultation=true`). The link cannot close while the
  segment between the two spacecraft meets the planet, a sphere of its
  equatorial radius. `occultation=(; min_ray_altitude_m=h)` raises that sphere
  by `h` metres, so that a beam grazing the atmosphere counts as blocked. The
  atmosphere's attenuation is not modelled.
- **Divergence** (`beam`). The force is scaled by
  `beam_force_factor(beam, link, terminal_a, terminal_b, d)` at separation `d`.
  - For a `GaussianBeam`, this is the share of the beam that the far aperture
    captures, `1 − exp(−2 r² / w(d)²)`, from the same beam functions as a
    power link (`beam_spot_radius_m`, `beam_capture_fraction`).
  - Within the Rayleigh range `π w0² / λ`, the force is close to the near-field
    force `magnification × power_w / c`. Far beyond it, the force falls as
    `1/d²`.
  - This is an engineering model, not a cavity model: it takes the
    magnification as the near-field value and scales it by one pass's
    geometric capture.
  - A beam model of your own subtypes `AbstractLaserBeam` and implements
    `beam_force_factor`, for example a recycling cavity or a measured curve.
- **Pointing** (`pointing=TerminalPointing(; field_of_regard_rad,
  error_rad=0.0)`).
  - The link can close only while each end's line of sight to the other lies
    within `field_of_regard_rad` of its terminal's boresight, read in the
    spacecraft's attitude. Within that cone the terminal points ideally: there
    is no slew rate or gimbal state.
  - A static pointing error `error_rad` scales the force by
    `exp(−2 (error / θ)²)`, with `θ = λ / (π w0)` the beam's divergence. It
    needs a `GaussianBeam`. Jitter is not modelled.
- **Mounted recoil** (`recoil=:mounted`). Each end's force acts at its
  terminal's mounted location (`terminal_pose`), so it turns the spacecraft with
  the torque of the terminal's lever arm, as the terminal's own recoil would.
  The default, `:central`, applies the force with no torque.

The solver does not locate occultation or field-of-regard boundary crossings
as events. Link availability and selection are recomputed at each accepted step;
a crossing between steps is not resolved at its exact time. Check timestep
convergence when these boundaries matter for your case.

Pointing and mounted recoil need terminal hardware at both ends (a
`LaserTerminal` carries a `LaserBeamLink`), and attitude (`orientation_sim=true`).
`register_link!` refuses a virtual terminal for them, and a run without attitude
refuses them at setup. With either option, range, occultation, divergence and
force direction all use the ray between the mounted terminals, including the
attitude transformation of each mount offset. Other configurations retain the
spacecraft body-reference positions. Call the terminal-aware `interlink_wrench`
and `link_available` forms for pointing or mounted recoil; their four-argument
force/availability forms lack the required terminal hardware.

The two forces are always equal and opposite, so the link conserves the run's
linear momentum. While the link is active, the terminals' own recoil applies
nothing: the link's force is their recoil (see
[Terminal hardware and power](#Terminal-hardware-and-power)). A checkpoint
records every option, and a resume with any of them changed is refused. The
gve example runs with a beam when given one:
`gve_sma_configuration(beam=GaussianBeam(...))`.

## Power links

A `PowerBeamLink` carries power, one way, from a `LaserTerminal` (the
transmitter) to a `PowerReceiver` (the receiver). It applies no force. This
configuration fragment uses the 50 W `terminal` mounted above and requests
40 W, within its limit:

```julia
receiver = PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4)
add_actuator!(other, other.root, receiver)
register_link!(model, spacecraft, (1, terminal), (2, receiver),
    PowerBeamLink(power_w=40.0, wavelength_m=1.064e-6, waist_radius_m=0.05, range_m=50e3))
args = SimulationConfiguration(...; interlink_model=model,
    scheduling_policy_model=PowerDeliveryPolicy())
```

- **The beam.** While the link is active, the transmitter runs at `power_w`
  optical watts. It draws `power_w / efficiency` from its own battery as a
  tier-1 laser load (see [Terminal hardware and power](#Terminal-hardware-and-power)).
  The beam is Gaussian:
  - its waist, of 1/e² radius `waist_radius_m`, is at the transmitter;
  - its spot radius at range `d`, `beam_spot_radius_m(link, d)`, is
    `w0 √(1 + (d / z_R)²)`, with the Rayleigh range `z_R = π w0² / λ`.
- **The receiver.** It captures `beam_capture_fraction(link, radius_m, d)`,
  `1 − exp(−2 r² / w²)`, of the beam through its aperture. It delivers that
  times its `efficiency` to its bus, up to its `max_power_w`
  (`link_received_power_w`).
  - For example, with a 1.064 µm beam from a 5 cm waist, a 25 cm receiver
    catches nearly all of it out to about 10 km, about a quarter at 100 km, and
    under 1 % at 1000 km.
- **Generation.** The delivered power is generation in the receiver's power
  budget. It shows in `power_generated_w` and in its own `power_received_w`
  column, charges the receiver's battery, and helps an empty battery recover.
  A full receiver battery curtails it, as it does solar power; the transmitter
  keeps paying.
- **Availability.** The link can close within `range_m`, while the segment
  between the spacecraft clears the planet (a sphere of its equatorial radius).
  When the transmitter's battery sheds tier 1, the beam becomes unavailable
  until the tier clears.
- **Direction.** A `PowerReceiver` can only receive, and a `LaserTerminal` can
  only transmit a power beam (`terminal_supports(terminal, link, side)`).
- **Scheduling.** Policies that score links by their forces never select a
  power beam, so use `PowerDeliveryPolicy()`, which scores each link by the
  power it would deliver, or `EnergyBalancingPolicy()`, which charges the
  lowest state of charge (see [Energy balancing](#Energy-balancing)). To weigh
  momentum and power links together, use a [weighted mix](#Weighted-mixes).
- **Columns.** A run with a power beam adds `power_received_w` and
  `power_beamed_w` per spacecraft to the [power columns](output_columns.md#Power-budget).
- **Not yet modelled.** Pointing and jitter, the receiver's incidence angle,
  the atmosphere, beam quality, and the beam's radiation pressure are not
  modelled.

## Data links

A `DataLink` carries data, one way, from a transmitter terminal to a receiver.
It applies no force. Its `budget` sets the rate:

```julia
radio() = RadioTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=5.0, gain_dbi=6.0, efficiency=0.5,
    noise_temperature_k=500.0)
sender, receiver = radio(), radio()
add_actuator!(sc, bus, sender)
add_actuator!(other, other.root, receiver)
link = DataLink(transmit_power_w=1.0, range_m=2000e3,
    budget=ShannonBudget(bandwidth_hz=1e6, frequency_hz=2.2e9, implementation_efficiency=0.5))
register_link!(model, spacecraft, (1, sender), (2, receiver), link)
args = SimulationConfiguration(...; interlink_model=model,
    scheduling_policy_model=DataThroughputPolicy())
```

- **The budget.**
  - `FixedRateBudget(; rate_bps)`: a constant rate at any range up to
    `range_m`. It works between virtual terminals too.
  - `ShannonBudget(; bandwidth_hz, frequency_hz, implementation_efficiency=0.5, max_rate_bps=Inf)`:
    a radio link budget between two `RadioTerminal`s.
    - The received power follows Friis: `P_r = P_t G_t G_r (λ / 4πd)²`, with
      the link's `transmit_power_w` as `P_t`, the two radios' `gain_dbi` and
      `λ = c / frequency_hz`.
    - The noise is `k_B T bandwidth_hz`, with the receiver's
      `noise_temperature_k` as `T`.
    - The rate is `implementation_efficiency × bandwidth_hz × log₂(1 + P_r / N)`,
      up to `max_rate_bps`.
    - For example, in the code above the rate is about 5.7 Mbit/s at 10 km,
      2.4 Mbit/s at 100 km and 0.17 Mbit/s at 1000 km.
  - A budget kind of your own subtypes `AbstractDataBudget` and implements
    `budget_rate_bps(budget, link, transmitter, receiver, distance_m)`. The
    terminals are `nothing` at a virtual end.
  - The rate is capped at each `RadioTerminal`'s `max_rate_bps`
    (`link_data_rate_bps`).
- **The radios.** A `RadioTerminal` can transmit or receive. While it transmits,
  it runs at the link's `transmit_power_w` and draws `transmit_power_w / efficiency`,
  as a tier-1 load. A receiving radio draws nothing; give its electronics a
  constant `PowerLoad`. When the transmitter's battery sheds tier 1, the link
  becomes unavailable until the tier clears. A `LaserTerminal` or a
  `PowerReceiver` cannot carry a data link.
- **Bits.** In a run with a data link, every spacecraft has two more integrated
  states, `data_sent_bits` and `data_received_bits`. Their rates are the sums of
  the rates its active links carry out and in, evaluated at every
  right-hand-side call. The solver integrates these counters alongside the
  other states, subject to its integration tolerances. They continue across a
  checkpoint resume.
- **Availability.** As for a power beam: within `range_m`, while the segment
  between the spacecraft clears the planet.
- **Direction.** Two-way traffic is two links, one each way. A radio carries
  one active link at a time, so the scheduler selects one direction of a
  two-way pair between the same two radios: the one with the higher score.
- **Scheduling.** Use `DataThroughputPolicy()`, which scores each link by the
  rate it would carry. The force and power policies never select a data link.
- **Columns.** A run with a data link adds `data_rate_tx_bps`,
  `data_rate_rx_bps`, `data_sent_bits` and `data_received_bits` per spacecraft
  (see [the data columns](output_columns.md#Data-links)).
- **Not yet modelled.** Pointing and antenna patterns, Doppler, atmospheric and
  polarization losses, interference between links and acquisition time are not
  modelled. Onboard storage and ground stations are below.

### Data storage

A spacecraft can hold data in a `DataStore` and produce it with data sources.
The data links then carry only what the stores hold, or what passes through
them:

```julia
add_data!(sc, bus, DataStore(capacity_bits=8e9, initial_bits=0.0))   # the root link, at most one
add_data!(sc, bus, DataSource(rate_bps=2e5))                          # any link, any number
```

- **States.** In a run with a store, every spacecraft has three more
  integrated states, which stay zero for a spacecraft without a store:
  - `data_stored_bits`, the data its store holds;
  - `data_produced_bits`, what its sources have produced;
  - `data_lost_bits`, the production a full store could not keep.

  They balance: `stored = initial + produced + received − sent − lost`. The
  counters `data_sent_bits` and `data_received_bits` count what the links
  actually carry.
- **Flow control.** A link carries at most its rate (`link_data_rate_bps`),
  and less at a store, so links lose nothing:
  - **An empty store** passes what comes in straight on: its production and
    what its incoming links carry, shared among its outgoing links in
    proportion to their rates. A relay with an empty store forwards data at
    zero stock.
  - **A full store** takes in no more than it sends on. Its production beyond
    that is lost and counted in `data_lost_bits`.
  - **A spacecraft without a store** neither limits its links nor is limited,
    as before.

  The rates are found at every right-hand-side evaluation, along chains of
  relays.
- **Full and empty.** A store counts as full within `1e-9` of its capacity and
  as empty within `1e-9` of its capacity of nothing. Events mark it full or
  empty when it crosses either, as a battery's do, and the flags are checked
  again at the end of each step. They continue across a checkpoint resume. A
  resume with another capacity or initial content is refused.
- **Sources.** `DataSource(; rate_bps, name="")` produces a constant rate while
  its spacecraft is active. A source kind of your own subtypes
  `AbstractDataSource` and implements
  `data_production_bps(source, x::StateSample, ctx::LinkContext, t)`, its rate
  in bit/s: for example imaging only over a region, a duty cycle, or a rate that
  depends on the attitude. A run refuses at setup a source on a spacecraft
  without a store, and a source kind without a method.
- **Solvers.** As with a battery, a run with a data store needs a solver mode
  that integrates the full right-hand side.
- **Scheduling by backlog.** `DataBacklogPolicy(; horizon_s)` scores a data
  link by what it could move in the next `horizon_s` seconds,
  `min(c, stored_tx / horizon_s + production_tx, room_rx / horizon_s + outflow_rx)` in bit/s,
  with `c` the link's rate and `production_tx` the transmitter's current source
  production rate. `outflow_rx` is the receiver's currently carried outgoing
  rate, sampled once from previously active links that remain available at the
  scheduling step, before choosing the new matching. It uses the runtime flow
  caps; it does not assume links in the new matching will be active. This is a
  scoring estimate only: actual transfers still obey the store's flow limits.
  It samples production at the current state and time and
  estimates that rate continuing over the horizon; it does not predict future
  source changes or incoming relay traffic. An empty sender that is still
  producing can keep its link active. A sender with neither stored data nor
  production, or a receiver with neither free space nor outflow, scores zero. A transmitter
  without a store counts as always holding data, and a receiver without one as
  always having room. It shares `DataThroughputPolicy`'s unit, so the two mix in
  a `WeightedPolicy` without rescaling. A policy of your own reads the stores
  with `stored_bits(ctx, i)` (NaN without a store) and `store_room_bits(ctx, i)`
  (Inf without one).
- **Backlog scoring limits in relay networks.** The score estimates a local
  transfer, not additional data delivered to its final destination. A full relay
  can attract incoming traffic that displaces its own production while a better
  direct path goes unused. Two full stores can also keep exchanging data without
  increasing delivery to a destination. Each incoming candidate independently
  receives the receiver's whole carried-outflow credit; that credit is neither
  reduced by the receiver's own production nor divided among candidates. If
  several links are selected, actual transfers still share the available flow
  capacity and respect store limits. These scoring limitations remain in v0.7.1;
  compare delivered and lost data when choosing a policy for a relay network.
- **Columns.** A run with a store adds `data_stored_bits`, `data_store_fill`
  (the fraction of the capacity held), `data_produced_bits` and
  `data_lost_bits` per spacecraft (see
  [the storage columns](output_columns.md#Data-storage)). `data_rate_tx_bps` and
  `data_rate_rx_bps` show the rates the links carry.
- **Not yet modelled.** One store holds one stream of data: there are no
  priority classes, no queues per destination, and no age or latency of the
  data held.

### Ground stations

A ground station is a fixed point on the rotating planet. Its antennas are link
terminals:

```julia
dish = RadioTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=20.0, gain_dbi=30.0, efficiency=0.5,
    noise_temperature_k=150.0)
station = GroundStation(name="north", lat_deg=42.3, lon_deg=-83.7, alt_m=270.0, min_elevation_deg=10.0,
    antennas=[dish])
s = add_ground_station!(model, station)                      # 1 for the first station
register_link!(model, spacecraft, (1, sender), (-s, 1), link)  # a downlink
# or a link from every spacecraft to every station:
register_links!(model, spacecraft, station_pairs(spacecraft, model), link)
```

- **Nodes.** Station `s` is the node `-s`, so its antenna `k` is the endpoint
  `(-s, k)` in link keys, the history, the contact windows and the link
  schedule. A link joins a spacecraft and a station; two stations cannot link.
  `forbidden_pairs` may name `(i, -s)`.
- **Where it is.** The station sits at its geodetic latitude, longitude and
  altitude on the planet's ellipsoid, and turns with the planet. Its inertial
  position is `L_PIᵀ r` and its velocity `L_PIᵀ (ω × r)`, with `r` its
  planet-fixed position, `L_PI` the planet's frame and `ω` its spin.
  - Availability and the scheduling policies read the frame at each accepted
    step.
  - The rates read it at the time of each right-hand-side evaluation, so a
    station moves through a step instead of holding its place at the step's
    start.
- **Availability.** A link to a station is available while the spacecraft is at
  least `min_elevation_deg` above the station's local horizon, the plane normal
  to its geodetic up direction, and within the link's `range_m`. The mask is the
  line of sight: a `DataLink` does not test the segment to a station against
  the planet.
- **Antennas.** They are `RadioTerminal`s, numbered `1:length(antennas)`, and
  carry one active link at a time. They steer freely above the mask, so their
  `boresight` is not read. A station has no power budget: its antennas draw
  nothing and are never shed.
- **Data.** A downlink counts in its spacecraft's `data_sent_bits`. An uplink, a
  link from `(-s, k)` to a spacecraft's radio, counts in that spacecraft's
  `data_received_bits`. In a run with data links, each station also counts the
  bits it has received and sent, two integrated states (`u.gs[s].received_bits`
  and `u.gs[s].sent_bits`) written as the `gs{s}_received_bits` and
  `gs{s}_sent_bits` columns (see [the data columns](output_columns.md#Data-links)).
  They continue across a checkpoint resume.
- **With data stores.** A station is a store with unbounded capacity and supply,
  so it never limits a link: a downlink carries what its spacecraft's store
  sends, and an uplink what its spacecraft's store can take in.
  `stored_bits(ctx, -s)` and `store_room_bits(ctx, -s)` are Inf, so
  `DataBacklogPolicy` drains the fullest stores to the ground first.
- **Scheduling.** Station links are candidates like any other. A station is never
  an orbital policy's target. It has no battery, so its `state_of_charge` is NaN
  and `EnergyBalancingPolicy` scores a link to it 0.
- **Outputs.** In the link schedule, a station end shows as `-s` in
  `a_spacecraft` or `b_spacecraft`, and by its registered name in `a_name` or
  `b_name`. These name columns are present only for station runs. Saved scenes
  include each site's planet-fixed coordinates; the station-capable Visualizer
  draws green markers and names on the rotating globe, and active or contact
  links in either direction. The labels toggle controls names and the interlinks
  toggle controls both stations and their links. Existing spacecraft-only scenes
  and schedule tables retain their format.
- **Checkpoints.** A run with stations records each station's site, mask and
  antennas, and refuses to resume when any of them has changed.
- **Not yet modelled.** Terrain and horizon masks, antenna steering limits and
  patterns, and atmospheric losses.

## Scheduling

At the start of the run and after every accepted solver step, the engine marks
a candidate available when its spacecraft pair is not forbidden, both
spacecraft are active, neither terminal is shed, and its link's
`link_available` allows it. The run's scheduling policy then scores each
available link, and the scheduler selects the terminal matching with the largest
total score. Only links with a positive score can be selected. The selected schedule is held during trial
stages, which avoids changing it during rejected steps. An inactive endpoint
contributes no force even if an earlier schedule still names its link.

In a run with interlinks, every spacecraft has two more integrated states.
`interlink_dv` integrates the acceleration its selected links' forces give
(m/s, inertial). `interlink_delta_sma` integrates their contribution to the
instantaneous osculating semimajor-axis rate (m); it is not a comparison against
a separately propagated unforced trajectory. Both stay zero for links that apply
no force, such as power beams and data links. They are saved in the
[interlink columns](output_columns.md#Interlinks). The history belongs to the
returned run configuration. With the default state isolation, the input
configuration is unchanged.

### Policies

A scheduling policy is a value of a policy type. Give the run one with
`SimulationConfiguration(scheduling_policy_model=policy)`; the default is
`SemimajorAxisRatePolicy()`. The run keeps it in a `SchedulingPolicyModel`, the
scheduler, which also holds the matching workspace:
`args.scheduling_policy_model.objective` is the policy.

| Policy | Scores each available link by | Unit |
| --- | --- | --- |
| `SemimajorAxisRatePolicy(; targets=1)` | the rate at which its force raises the osculating semimajor axis of its target spacecraft, `2a²(v·F)/(μm)`, summed over the target ends | m/s |
| `EccentricityRatePolicy(; targets=1)` | the rate at which its force raises its target's eccentricity, `ê·ė` | 1/s |
| `PowerDeliveryPolicy()` | the electrical power it would deliver to the receiver's bus, `link_received_power_w` | W |
| `DataThroughputPolicy()` | the data rate it would carry, `link_data_rate_bps` | bit/s |
| `DataBacklogPolicy(; horizon_s)` | the data it could move in the next `horizon_s`: the least of its rate, the transmitter's stored data over `horizon_s` plus its current production rate, and the receiver's free space over `horizon_s` plus current carried outflow | bit/s |
| `EnergyBalancingPolicy()` | the power it would deliver times the state-of-charge gap from sender to receiver, `P × (s_tx − s_rx)` | W |
| `WeightedPolicy(terms...)` | `Σ weight × (score / scale)` over its terms | 1 |

The orbital policies take each link's force from its `link_force`. For
`SemimajorAxisRatePolicy`, `targets` is a spacecraft vector index, a vector or
range of them, `:both` (both ends of every link) or `:lower_index` (the lower
spacecraft index of each pair). For `EccentricityRatePolicy` it is an index or
`:both`. A link that lowers the objective scores negative, so it is never
selected.

### Energy balancing

`EnergyBalancingPolicy()` charges the lowest state of charge. It scores a power
beam into terminal hardware, between two spacecraft with a battery, by the
power it would deliver times the gap between the two batteries' states of
charge, `P × (s_tx − s_rx)`. So:

- power flows only from a fuller battery to an emptier one, and a beam into a
  full battery never scores: a battery the power system counts as full (its
  stored energy within `1e-9` of its capacity) counts as 1, so two full
  batteries tie at 0;
- of receivers competing for one sender's terminal, the emptier one wins at
  equal delivered power, while a much better-coupled receiver can still win
  over a slightly emptier one;
- the score is signed: a beam from an emptier battery to a fuller one scores
  negative, so it is never selected alone and counts against the beam in a mix;
- a link with no battery at either end, or one that would deliver no power,
  scores 0.

The selection changes only at accepted steps, so when two states of charge
cross within a step, the beam runs to the step's end. `PowerDeliveryPolicy`,
by contrast, beams into a full battery: the receiver curtails, and the sender
still pays.

### Weighted mixes

`WeightedPolicy(PolicyTerm(policy; scale, weight=1.0)...)` scores each link by
`Σ weight × (score / scale)` over its terms, summed in term order. `scale` is in
the term's own unit, the amount of its score that counts as one, so each term is
dimensionless before the weights apply; `weight` is finite and may be negative.
The scales are fixed, not normalized by each step's best link, so a link's
score does not depend on the other candidates.

```julia
policy = WeightedPolicy(
    PolicyTerm(SemimajorAxisRatePolicy(targets=1); scale=1.0),   # 1 m/s of SMA rate counts as one
    PolicyTerm(EnergyBalancingPolicy(); scale=100.0))            # and so do 100 W
args = SimulationConfiguration(...; interlink_model=model, scheduling_policy_model=policy)
```

Each term applies its own rules: an orbital term scores power beams 0 and a
power term scores momentum links 0, so a mix weighs link kinds against each
other. Each term's checks run, and an eccentricity term keeps its near-circular
refusal over the mix's candidates. Mixes nest.

### A custom policy

A policy type subtypes `AbstractSchedulingPolicy` and implements
`candidate_score`. This one prefers short links:

```julia
struct ShortLinkPolicy <: SpaceAGORA.AbstractSchedulingPolicy end
SpaceAGORA.candidate_score(::ShortLinkPolicy, candidate, ctx, t) =
    1 / sqrt(sum(abs2, candidate.a.pos_ii - candidate.b.pos_ii))

args = SimulationConfiguration(...; interlink_model=model, scheduling_policy_model=ShortLinkPolicy())
```

- **The candidate** (`LinkCandidate`) has its `key`, its `link`, its endpoints'
  state samples `a` and `b` (in key order, so a directed link's transmitter is
  `a`), and the terminal hardware at its ends, `terminal_a` and `terminal_b`
  (`nothing` for a virtual terminal).
- **The context** (`SchedulingContext`) has the planet's gravitational parameter
  `mu`, the `link_context` to pass to the link hooks, and the run's `spacecraft`;
  `state_of_charge(ctx, i)` gives spacecraft `i`'s battery state of charge (NaN
  without a battery), and `stored_bits(ctx, i)` and `store_room_bits(ctx, i)`
  the data its store holds and the room left (NaN and Inf without a store). `t`
  is the elapsed time in seconds.
- **The rules.** The hook is called for available candidates only; the others
  score 0. A nonfinite score is refused, and only positive scores can be
  selected. A run refuses at setup a policy type without a `candidate_score`
  method.
- **Plain values.** A policy holds no run state: it is a value whose fields say
  what it scores. A checkpoint records its type and every field by value
  (numbers with their bits; arrays, tuples, named tuples and structs element by
  element and field by field), so an equal policy rebuilt for the resume
  resumes, and a resume with another value is refused. A function field is
  recorded by its type and, for a closure or a callable struct, every value it
  holds; a named function is recorded by its type alone, so a change to its
  code is not detected. A closure reads the globals it uses when it is called,
  and their values are not recorded: build it in a function or a `let` block,
  so that it captures them. A checkpointed run refuses at setup what a
  checkpoint cannot record by value: a closure that captures a variable
  reassigned after the closure was made, a dictionary or a set, a stream, a
  task or a lock, a pointer or a value holding one (such as a `Regex`), whose
  address differs from one process to the next, an unassigned field or element,
  and fields that contain themselves.

### Contact windows

The history (`model.history`) has one `InterLinkHistory` sample per accepted
step: its time, the links the scheduler selected (`active`) and the links that
were available (`available`). The available links are recorded in a run with a
data link, or when you ask for them with `InterlinkModel(; record_available=true)`;
`record_available=false` leaves them out. In a large constellation over a long
run they are a list per step, so other runs leave them out by default. `contact_windows(history)` turns it into windows,
one `(key, start_s, end_s, duration_s)` for each span over which a link could
close; `contact_windows(history; active=true)` gives the spans over which it
was selected. `contact_time_s(history, key)` sums a link's windows.

```julia
solution = run_simulation(args; return_solution=true)
history = solution.prob.p.args.interlink_model.history
windows = contact_windows(history)
contact_time_s(history, ((1, 1), (2, 1)); active=true)
```

A window starts at the first accepted step that has the link and ends at the
next step that does not, or at the last step: the resolution is the
scheduler's. A checkpointed run keeps its history in
`simulation_checkpoint.interlinks.bin` beside its checkpoint, and a
resumed run takes it back, so a resumed history holds the whole run, as the
results do (see [Outputs](outputs.md)). If that record is missing or unreadable,
the run warns and records only the resumed interval. A resume that solves nothing
returns an empty schedule in that case. With the same configuration, checkpoint
schedule and deterministic models whose complete run state is restored, resumed
history matches the uninterrupted checkpointed run. A run without checkpoints
can take different adaptive steps and record different scheduling samples.

### The link schedule

Each history sample also records the score that selected each active link
(`scores`), read from the scheduler after it selects. `link_schedule` turns the
history into a table with one row per window: a span over which a link was
selected (`window == "active"`) and, when the history records the available
links, a span over which a link could close (`"contact"`). Each row gives the
link (its number, kind, direction and endpoints), the span (`start_s`, `end_s`,
`duration_s` and the number of samples, `steps`) and, for an active row, the score
at its first and last sample and its extremes.

```julia
results = run_simulation(args; return_results=true)
schedule = link_schedule(results)              # or link_schedule(model; policy=scheduler)
active = schedule[schedule.window .== "active", :]
```

The scores are in the policy's unit (see [Policies](#Policies)): m/s for
`SemimajorAxisRatePolicy`, W for the power policies, bit/s for
`DataThroughputPolicy`, 1 for a `WeightedPolicy`. The table's metadata names the
policy and the unit. With `InterlinkModel(; record_terms=true)` and a
`WeightedPolicy`, each sample also records each top-level term's own score, in
its own unit, and the table gains `term1_start`, `term2_start`, ...: each term's
score at the span's first sample. A run with another policy refuses
`record_terms=true`. Recording the terms leaves the mix's scores, and so the
schedule, unchanged.

A run that writes results writes the table beside them, as
`simulation_results_link_schedule.csv` (when `save_csv` is set) and
`simulation_results_link_schedule.feather`, and lists both in the results
manifest with the policy and its unit (see the
[link schedule's columns](output_columns.md#Link-schedule)). The table is not
thinned with the results rows. A run without interlinks writes none, and
`link_schedule(results)` returns the table's columns with no rows. Like the
contact windows, the schedule has the scheduler's resolution. A resumed run
restores its recorded history from beside the checkpoint, including earlier
windows. Its metadata and manifest give `history_start_s` and `history_end_s`.
If the history is unavailable, a warning identifies the missing earlier interval;
that interval means not recorded, rather than no link. With `save_visualization_scene`, the viewer payload
carries the schedule's windows, so the viewer can draw each link while it is
active (see [the viewer page](outputs.md#Interactive-viewer-page)).

### Optional matching backend

Tree-shaped components are solved in Julia. Of parallel links (directed links
both ways between the same two terminals), only the higher-scoring one is a
candidate, so a two-way pair alone is a tree. Cyclic components use NetworkX's
maximum-weight matching through the optional PythonCall extension. Core loading
and tree matching do not initialize Python. A cyclic component without the
backend produces an installation error rather than using an approximate match.

Install PythonCall in the Julia environment used for the run, and install
NetworkX (`>=3.4,<4`) in a compatible Python environment. Set
`JULIA_PYTHONCALL_EXE` to that interpreter before starting Julia; when managing
Python yourself, also set `JULIA_CONDAPKG_BACKEND=Null`. Then explicitly load
`PythonCall` and `SpaceAGORA` in either order. For example, in a separate Julia
project containing SpaceAGORA:

```julia
using Pkg
Pkg.add("PythonCall")
using SpaceAGORA, PythonCall
```

The backend requires Python's global interpreter lock. Run cyclic matching on
its owning thread, or use independent process workers. Merely adding more Julia
threads does not make arbitrary calls into Python safe. Both matching paths use
exact integer weights derived from the floating-point scores. Equal-score ties
may choose different valid optima across backend versions.

## A custom link kind

A kind implements `link_range_m` and the effects it has. This one pulls its
endpoints together with a constant force:

```julia
struct PullLink <: SpaceAGORA.AbstractInterlink
    force_n::Float64
    range_m::Float64
end
SpaceAGORA.link_range_m(link::PullLink) = link.range_m
function SpaceAGORA.link_force(link::PullLink, a, b, ctx)
    direction = (b.pos_ii - a.pos_ii) / sqrt(sum(abs2, b.pos_ii - a.pos_ii))
    return (link.force_n * direction, -link.force_n * direction)
end

register_link!(model, spacecraft, (1, 1), (2, 1), PullLink(3.0, 2_000.0))
```

A run refuses at setup a link whose kind has no `link_range_m` method. Links of
different kinds can be registered in one model.

At a ground station end (see [Ground stations](#Ground-stations)), the hooks
receive the station's state sample: the station is its `spacecraft`, its
`spacecraft_index` is `-s`, and it has no mass and no attitude. The station sits
on the planet, so a kind that tests the planet's occultation should leave a
station end to the station's elevation mask, as `DataLink` does
(`a.spacecraft isa GroundStation`). A station's antennas are `RadioTerminal`s,
which carry a kind only where `terminal_supports` allows it.

### Terminals and torques

A kind that reads its terminal hardware, or that turns the spacecraft it
links, implements the terminal-aware forms of the hooks. The engine calls these
forms, and their defaults call the ones above, so a kind written without them
works as before, bit for bit.

- `link_available(link, terminal_a, terminal_b, a, b, ctx)` decides whether the
  link can close, given the terminal hardware at its two ends (`nothing` for a
  virtual terminal).
- `interlink_wrench(link, terminal_a, terminal_b, a, b, ctx)` returns
  `(force_a_ii, torque_a_body, force_b_ii, torque_b_body)`.
  - The forces act as `link_force`'s do, and the scheduling policies score
    candidate links by them.
  - Each torque is in its spacecraft's body frame, about its body reference
    point, as a mounted actuator's torque is (see
    [Frames](../reference/extension_conventions.md#Frames)).
  - The engine adds each torque to its spacecraft's angular acceleration as
    `I⁻¹ τ`, with `I` the inertia tensor. This is exact for the rigid spacecraft
    interlinks support.
- `link_applies_torque(link) = true` says that a kind's torques can be nonzero.
  A run with such a link needs attitude (`orientation_sim=true`), and setup
  refuses it otherwise. A nonzero torque from a kind that does not say so
  stops a run without attitude.
- `terminal_pose(spacecraft, terminal)` gives a mounted terminal's boresight and
  location in its spacecraft's body frame. A force applied at that location has
  the torque the terminal's own recoil would have: the cross product of the
  location and the force, in the body frame.

### Supported scope

Interlinks apply forces, and through `interlink_wrench` torques, to rigid
spacecraft, or carry power or data. The laser momentum link does not model terminal pointing,
occultation, beam losses or torques; `LaserBeamLink` adds them as options (see
[Laser physics](#Laser-physics)). Neither models the terminals' thermal
dynamics. Links through hardware terminals draw battery
power and close when it is shed; links between virtual terminals draw none.
Actuator faults on a link-driven terminal do not yet act on its link. The
spacecraft's `battery_energy_index` and `tempurature_index` fields are
metadata, not availability constraints. The cavity
magnification and force law are modeling assumptions; passing implementation
tests does not validate a physical laser system.

Use a full-state solver such as `:tsit5`, or the opt-in `:multirate_mri` route
described in [Solver configuration](solver_configuration.md#Coupled-multirate-integration).
The original split modes (`:split_imex`, `:multirate`, `:gravity_backbone_split`),
external propagator coupling, articulated/compliant vehicles and independent
constellation ensembles are rejected for linked spacecraft. Near circularity, `EccentricityRatePolicy`
also rejects candidate sets that can simultaneously apply several positive-score
links to one target: individual eccentricity-rate norms are not additive there.
`SemimajorAxisRatePolicy` does not have that restriction.

Scheduling follows accepted steps, so changes of range or preferred link are
resolved at that cadence. Check timestep convergence for a new case. A checkpoint
resume must retain the interlink topology, link kinds (each named by its module
and type, with its concrete type parameters) and parameters, scheduling policy
(its type and every field), and each spacecraft's terminal hardware: its kinds,
in mounting order, and every field of each, such as a receiver's aperture or a
transmitter's efficiency; a new run is required when those change. Checkpoints
written by SpaceAGORA 0.4 and earlier still resume, and so do checkpoints written
with the symbol policies of 0.5, under the typed policy that schedules the same.
A resumed run takes its interlink history back from beside its checkpoint, as
its results take back their rows, so its history, schedule and results table
hold the whole run, and integrated diagnostic states continue from the
checkpoint. Restoring recorded history requires the same `record_available` and
`record_terms` settings; changing either is refused before result streams are
changed. Start a new run to change these settings.

## Deprecated names

The names from when every link was a laser link keep working until SpaceAGORA
1.0. The old keywords, fields and functions warn once per process, and throw
under `julia --depwarn=error`. The two old type names are aliases of the new
types and cannot warn on their own: `InterLinkParameters()`,
`InterLinkParameters(terminal)`, the positional `InterLinkParameters(P, B, range)`
and `ConstellationInterLinkModel(; forbidden_pairs)` run silently, so search for
those names when you move to the new ones. The `laser_dv` and `laser_delta_sma`
columns cannot warn either: a run writes them beside the `interlink_dv` and
`interlink_delta_sma` columns, with the same values, until 1.0. Leave the
`:laser_dv` and `:laser_delta_sma` fields out of an explicit `save_fields` list
to stop writing them.

| Deprecated | Use instead | Warns |
| --- | --- | --- |
| the columns `sc1_laser_dv_1` .. `sc1_laser_dv_3` and `sc1_laser_delta_sma`, and their save fields `:laser_dv` and `:laser_delta_sma` | `sc1_interlink_dv_1` .. `sc1_interlink_dv_3` and `sc1_interlink_delta_sma`; `:interlink_dv` and `:interlink_delta_sma` | no, a column |
| `ConstellationInterLinkModel` | `InterlinkModel` | no, an alias |
| its `constellation_link_type` keyword and property | nothing: each link has its own kind | yes |
| `InterLinkParameters` | `LaserMomentumLink` | no, an alias |
| its keywords and fields `P`, `B`, `range` | `power_w`, `magnification`, `range_m` | yes |
| `register_allowed_link!(model, spacecraft, a, b; parameters)` | `register_link!(model, spacecraft, a, b, link)` | yes |
| `force_on_endpoint(current, other, link)` | `link_force(link, a, b, ctx)` | yes |
| `laser_force_on_spacecraft(model, state, satellite)` | `interlink_force_on_spacecraft(model, state, satellite)` | yes |

## Removed in 0.6.0

The scheduling policies chosen by a symbol were removed in SpaceAGORA 0.6.0,
without a deprecation period, with their `target_index` and `target_indices`
keywords and properties and the older `target_idx` and `target_idxs`. A call
with a symbol throws an `ArgumentError` that names its typed replacement, also
under `julia --depwarn=error`, and so does reading or setting a removed
property, naming `.objective` or `.objective.targets`. Checkpoints written with
them still resume under the typed policy that schedules the same.

| Removed | Use instead |
| --- | --- |
| `SchedulingPolicyModel(:gve_sma)` | `SemimajorAxisRatePolicy()` |
| `SchedulingPolicyModel(:gve_sma; target_index=k)` | `SemimajorAxisRatePolicy(targets=k)` |
| `SchedulingPolicyModel(:gve_sma; target_index=nothing)` | `SemimajorAxisRatePolicy(targets=:both)` |
| `SchedulingPolicyModel(:gve_sma; target_indices=v)` | `SemimajorAxisRatePolicy(targets=v)` (a `target_index` beside it was ignored) |
| `SchedulingPolicyModel(:gve_sma_first_index)` | `SemimajorAxisRatePolicy(targets=:lower_index)` |
| `SchedulingPolicyModel(:gve_eccentricity; target_index=k)`, or `=nothing` | `EccentricityRatePolicy(targets=k)`, or `targets=:both` |
| `SchedulingPolicyModel(:power_delivery)` | `PowerDeliveryPolicy()` |
| `SchedulingPolicyModel(:data_throughput)` | `DataThroughputPolicy()` |
| the keywords `target_idx`, `target_idxs` | as `target_index`, `target_indices` above |
| `scheduler.policy` (a symbol) | `scheduler.objective` (a policy value) |
| `scheduler.target_index`, `scheduler.target_indices`, `scheduler.target_idx`, `scheduler.target_idxs` | `scheduler.objective.targets` |

The two integrated states the interlinks add to each spacecraft were renamed in
0.6.0 as well. The name of a state component cannot warn, so they were renamed
without a deprecation period; their columns keep the old names as aliases
until 1.0 (see [Deprecated names](#Deprecated-names)). Checkpoints keep the old
names inside, so checkpoints written before 0.6.0 resume.

| Removed | Use instead |
| --- | --- |
| the state `laser_dv`, as in `solution.u[end].sc[i].laser_dv` | `interlink_dv` |
| the state `laser_delta_sma` | `interlink_delta_sma` |

Checkpointed pointing and mounted-recoil beams validate the terminal's mounting
link as well as its hardware: changing its position, fixed attitude or hinge
axis/limits refuses resume. A `PanelHinge` restores its moving angle from the
checkpoint. A terminal on a panel driven by `SolarPanelAngleOfAttackControlModel`
needs a `PanelHinge` for checkpointing; arbitrary custom changes to a mount's
position or attitude need their own checkpoint state restoration. The root
link's attitude remains part of the integrated spacecraft state.

The schedule is also saved for long-format and checkpoint-streamed results. A
resume that only finishes interrupted trajectory file writes can rebuild its
schedule from validated checkpoint history. Without that history, it preserves
any earlier schedule file but does not list it in the recovered manifest or
returned files. A completed resume that solves nothing returns validated history
without rewriting files. An unusable history sidecar is preserved under a
`before_resume_at_<t>s` name immediately before the first replacement checkpoint;
a resume that fails before then leaves it untouched.
