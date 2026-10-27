# Interlinks

This page covers connecting spacecraft through interlinks and choosing, with a
scheduling policy, which of the links are active at each step.

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [Parallel Execution](parallel_execution.md)

## Interlinks and constellations

An interlink is a connection between two terminals that acts while the geometry
allows it: it pushes the two spacecraft apart, carries power, or carries data.
A constellation is a set of spacecraft and their orbits, laid out and run
together. The two are independent. A constellation needs no interlinks, and two
spacecraft can be linked without forming one.

## Links and link kinds

Each candidate link has a kind, which is a value of a subtype of
`AbstractInterlink`. The kind decides when the link can close and what it does
to its two endpoints, through three hooks:

- `link_range_m(link)`: the maximum separation in meters. Every kind
  implements it.
- `link_available(link, a, b, ctx)`: whether the link can close. The default is
  true while the endpoints are apart and within `link_range_m`.
- `link_force(link, a, b, ctx)`: the inertial forces (N) on `a` and `b` while
  the link is active. The default is no force.

`a` and `b` are the endpoints' state samples (`StateSample`, with position,
velocity, mass, attitude and `spacecraft_index`), in the order of the link's
key: the lower `(spacecraft, terminal)` pair first or, for a directed link, the
transmitter first. `ctx` is a `LinkContext` holding the time (s into the run)
and the run's planet.

The built-in kinds are:

- `LaserMomentumLink(; power_w=10_000.0, magnification=100.0, range_m=200e3)`:
  an open-cavity laser link with the idealized central force
  `magnification * power_w / c` on each spacecraft, directed away from the
  other, equal and opposite. `LaserMomentumLink(terminal)` takes its numbers
  from a `LaserTerminal`.
- `LaserBeamLink(; power_w, magnification, range_m, ...)`: the same link with
  optional occultation, beam divergence, pointing and mounted recoil (see
  [Laser physics](#Laser-physics)).
- `PowerBeamLink(; power_w, wavelength_m, waist_radius_m, range_m)`: a power
  beam from a `LaserTerminal` to a `PowerReceiver` (see
  [Power links](#Power-links)).
- `DataLink(; transmit_power_w, budget, range_m)`: a one-way data link from a
  `RadioTerminal` or a virtual terminal to another, at the rate its link budget
  sets (see [Data links](#Data-links)).

## Setting up interlinks

An `InterlinkModel` holds the candidate links, the spacecraft pairs that may
never link (`forbidden_pairs`), and the history of the links the scheduler
selected. An endpoint is a `(spacecraft vector index, terminal number)` pair,
where the index is the spacecraft's position in `dynamics_model.spacecraft`,
not its `id`. Each terminal takes part in at most one selected link at a time.

In this fragment spacecraft 1 gets a candidate link to each of spacecraft 2
and 3; `spacecraft` is the run's spacecraft vector and `...` stands for the
rest of the configuration:

```julia
model = InterlinkModel()
for partner in 2:3
    register_link!(model, spacecraft, (1, 1), (partner, 1),
        LaserMomentumLink(power_w=10_000.0, magnification=100.0, range_m=200e3))
end
args = SimulationConfiguration(...; interlink_model=model,
    scheduling_policy_model=SemimajorAxisRatePolicy())
```

The runnable starting point is `examples/gve_sma_interlinks.jl`. It builds
three spacecraft, registers two candidate links from spacecraft 1, and at each
step selects the one that raises spacecraft 1's semimajor axis faster.
`examples/targeted_40sat_sma_interlinks.jl` uses the same interface with thirty
helper spacecraft and ten targets. Both use inverse-square gravity and no
atmosphere, so they need neither native GRAM nor private telemetry. From the
repository root:

```julia
include("examples/gve_sma_interlinks.jl")
args = gve_sma_configuration(duration=120.0)
solution = run_simulation(args; return_solution=true)
history = solution.prob.p.args.interlink_model.history
```

## Terminals

A spacecraft's terminals are numbered in two groups: first its `n_terminal`
virtual terminals, `1:n_terminal`, then the terminal hardware mounted on it, in
mounting order. `n_terminal` defaults to 1, and a spacecraft with only hardware
terminals sets `n_terminal=0`. A `LaserTerminal`, `PowerReceiver` or
`RadioTerminal` mounted with `add_actuator!` is terminal hardware; a kind of
your own opts in with `SpaceAGORA.link_terminal(::MyTerminal) = true`.
`terminal_count(spacecraft)` gives the total, and
`mounted_link_terminals(spacecraft)` lists the hardware with its numbers.
Virtual terminals remain useful for idealized studies: they draw no power and
are never shed.

An endpoint's terminal can be given as a number or as the mounted terminal
itself. In this fragment the second spacecraft, `other`, is built the same way
as `sc`:

```julia
bus = Link(root=true, m=100.0)
sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, n_terminal=0)
terminal = LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=50.0, magnification=100.0, range_m=200e3)
add_actuator!(sc, bus, terminal)
# ... a second spacecraft, `other`, likewise
spacecraft = [sc, other]
register_link!(model, spacecraft, (1, terminal), (2, 1), LaserMomentumLink(terminal))
```

A terminal must be able to carry its link (`terminal_supports(terminal,
link)`), and `register_link!` refuses one that cannot. A `LaserTerminal`
carries `LaserMomentumLink` and `LaserBeamLink` and transmits a
`PowerBeamLink`. A virtual terminal carries any kind except one that reads the
hardware at both ends: a `DataLink` with a `ShannonBudget`, or a
`LaserBeamLink` with pointing or mounted recoil. `terminal_number(spacecraft,
terminal)` gives a mounted terminal's number. Two mounted terminals with equal
fields cannot be told apart, so for those give the number.

### Directed links

A kind with `link_directed(link) = true` runs one way, from the first endpoint
given to `register_link!` (the transmitter) to the second (the receiver). Its
key keeps that order, so its hooks receive the transmitter as `a`. The opposite
direction between the same terminals is a separate link. Power beams and data
links are directed; laser momentum and laser beam links are not.

### Registering many links

`register_links!(model, spacecraft, pairs, link)` registers one link per
spacecraft pair. `all_pairs(spacecraft)` gives every pair, or you can pass any
vector of `(i, j)` index pairs. For a constellation laid out by
`build_constellation` or `walker_constellation` (see
[Constellations](constellations.md)), the spacecraft's planes and slots give
the usual topologies:

- `same_plane_pairs(spacecraft)`: every pair within each plane.
- `ring_pairs(spacecraft)`: each slot with the next in its plane, closing the
  ring.
- `grid_pairs(spacecraft; wrap_planes=true)`: the "+grid", which is the rings
  plus each slot with the same slot of the next plane. Pass
  `wrap_planes=false` for a Walker star, whose first and last planes
  counter-rotate.

Spacecraft without a constellation slot take part in none of these pairs, and
two spacecraft holding the same slot are refused.

`link` is a kind value, or a function `(i, j) -> link` that builds each pair's
link. Each spacecraft's links take, in turn, the terminals that can carry them,
starting again from the first when they run out. Links that share a terminal
compete for it when the scheduler selects. Pass
`terminals=(i, j) -> (terminal_i, terminal_j)` to choose the terminals
yourself. For a directed link, `i` is the transmitter.

```julia
register_links!(model, spacecraft, all_pairs(spacecraft), LaserMomentumLink(power_w=40.0))
# or, for a constellation's +grid:
register_links!(model, spacecraft, grid_pairs(spacecraft), LaserMomentumLink(power_w=40.0))
```

### Terminal hardware and power

A hardware terminal that is a link endpoint is driven by its links for the
whole run. While one of its links is active, it runs at the level the link
sets, `link_terminal_levels(link)`; otherwise its level is 0. For a
`LaserMomentumLink` the level is the link's `power_w` at both ends. At that
level the terminal draws `actuator_consumption(terminal, level)` from its
spacecraft's power budget, which for a `LaserTerminal` is
`power_w / efficiency`. The draw appears in the power columns and in the
`actuator_level` and `actuator_power` columns.

While a link is active, only the link's force acts. The terminal's own recoil
applies only when a controller drives a terminal that is no link's endpoint. A
mounted-actuation controller's levels for a link-driven terminal apply nothing,
and setup logs a message saying so once. A link that would drive its terminal
past its limits (above a `LaserTerminal`'s `max_power_w`) is refused at setup;
build the link from the terminal, `LaserMomentumLink(terminal)`, to match
them.

When the battery sheds tier 1 (see [Powering a Spacecraft](power_budget.md)),
it sheds laser and radio terminals, and links through a shed terminal are
unavailable until the tier clears. A shed terminal keeps asking for the
highest level its links set, so the tier clears when the battery can supply
that power again, not merely because the link closed. A terminal kind's tier
is `terminal_shed_tier(terminal)`: 1 for `LaserTerminal` and `RadioTerminal`,
and 0 (never shed) for `PowerReceiver` and by default. Virtual terminals draw
nothing and are never shed.

## Laser physics

`LaserBeamLink` is the open-cavity laser link with four effects you can switch
on. With all four off it gives the same result as `LaserMomentumLink`, bit for
bit.

In this fragment `terminal_1` and `terminal_2` are `LaserTerminal`s mounted on
spacecraft 1 and 2, each with a `max_power_w` of at least 10 kW, and the run
has attitude:

```julia
beam = GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.05, aperture_radius_m=0.08)
link = LaserBeamLink(power_w=10_000.0, magnification=100.0, range_m=200e3;
    occultation=(; min_ray_altitude_m=100e3),   # or true
    beam,
    pointing=TerminalPointing(field_of_regard_rad=deg2rad(30)),
    recoil=:mounted)
register_link!(model, spacecraft, (1, terminal_1), (2, terminal_2), link)
```

- Line of sight (`occultation=true`). The link cannot close while the segment
  between the two spacecraft meets the planet, taken as a sphere of its
  equatorial radius. `occultation=(; min_ray_altitude_m=h)` raises that sphere
  by `h` meters, so that a beam grazing the atmosphere counts as blocked.
  Atmospheric attenuation is not modeled.
- Divergence (`beam`). The force is scaled by
  `beam_force_factor(beam, link, terminal_a, terminal_b, d)` at separation
  `d`. For a `GaussianBeam(; wavelength_m=1.064e-6, waist_radius_m,
  aperture_radius_m)` the factor is the share of the beam the far aperture
  captures, `1 − exp(−2 r² / w(d)²)`, from the same beam functions a power
  link uses (`beam_spot_radius_m`, `beam_capture_fraction`). The factor stays
  near 1 while the spot is smaller than the aperture, which for an aperture
  wider than the waist holds out to roughly the Rayleigh range `π w0² / λ`;
  there the force is close to the near-field force
  `magnification × power_w / c`. Far beyond it the force falls as `1/d²`. This
  is an engineering model, not a cavity model: it takes the magnification as
  the near-field value and scales it by one pass's geometric capture. A beam
  model of your own, such as a recycling cavity or a measured curve, subtypes
  `AbstractLaserBeam` and implements `beam_force_factor`.
- Pointing (`pointing=TerminalPointing(; field_of_regard_rad, error_rad=0.0)`).
  The link can close only while each end's line of sight to the other lies
  within `field_of_regard_rad` (at most π) of its terminal's boresight, in the
  spacecraft's attitude. Within that cone the terminal points ideally, with no
  slew rate or gimbal state. A static pointing error `error_rad` scales the
  force by `exp(−2 (error / θ)²)`, where `θ = λ / (π w0)` is the beam's
  divergence; it needs a `GaussianBeam`. Jitter is not modeled.
- Mounted recoil (`recoil=:mounted`). Each end's force acts at its terminal's
  mounted location (`terminal_pose`), so it also turns the spacecraft with the
  torque of the terminal's lever arm, as the terminal's own recoil would. The
  default, `recoil=:central`, applies the force with no torque.

The solver does not locate occultation or field-of-regard boundary crossings
as events. Availability and selection are recomputed at each accepted step, so
a crossing between steps takes effect at the next one. Check timestep
convergence when these boundaries matter for your case.

Pointing and mounted recoil need terminal hardware at both ends (a
`LaserTerminal` carries a `LaserBeamLink`) and attitude
(`orientation_sim=true`). `register_link!` refuses a virtual terminal for
them, and a run without attitude refuses them at setup. With either option,
range, occultation, divergence and force direction all use the ray between the
mounted terminals, with each mount offset rotated by its spacecraft's
attitude. Without them the link uses the spacecraft body-reference positions.
If you call the hooks yourself for such a link, use the terminal-aware forms,
`interlink_wrench` and the six-argument `link_available`. The four-argument
`link_force` and `link_available` have no terminal hardware to read, and
throw.

The two forces are always equal and opposite, so the link conserves the run's
linear momentum. While the link is active the terminals' own recoil applies
nothing: the link's force is their recoil (see
[Terminal hardware and power](#Terminal-hardware-and-power)).

A checkpoint records every option, and a resume with any of them changed is
refused. For pointing and mounted recoil the checkpoint also records each
terminal's mounting link, and a change to its position, its fixed attitude, or
its hinge axis or limits refuses the resume. A `PanelHinge` restores its angle
from the checkpoint. A terminal on a panel that
`SolarPanelAngleOfAttackControlModel` drives needs a `PanelHinge` in a
checkpointed run. Any other way of moving a mount during a run needs its own
checkpoint state restoration. The root link's attitude is part of the
integrated spacecraft state, so it needs nothing extra.

`gve_sma_configuration(beam=GaussianBeam(...))` runs the gve example with laser
beam links, which then also have `occultation=true`.

## Power links

A `PowerBeamLink` carries power one way, from a `LaserTerminal` (the
transmitter) to a `PowerReceiver` (the receiver), and applies no force. This
configuration fragment uses the 50 W `terminal` mounted above and asks for
40 W, within its limit:

```julia
receiver = PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4)
add_actuator!(other, other.root, receiver)
register_link!(model, spacecraft, (1, terminal), (2, receiver),
    PowerBeamLink(power_w=40.0, wavelength_m=1.064e-6, waist_radius_m=0.05, range_m=50e3))
args = SimulationConfiguration(...; interlink_model=model,
    scheduling_policy_model=PowerDeliveryPolicy())
```

While the link is active, the transmitter runs at `power_w` optical watts and
draws `power_w / efficiency` from its own battery as a tier-1 laser load (see
[Terminal hardware and power](#Terminal-hardware-and-power)). The beam is
Gaussian, with its waist, of 1/e² radius `waist_radius_m`, at the transmitter.
Its spot radius at range `d` is `beam_spot_radius_m(link, d)`,
`w0 √(1 + (d / z_R)²)`, with the Rayleigh range `z_R = π w0² / λ`.
`wavelength_m` defaults to 1.064 µm and `waist_radius_m` to 5 cm.

The receiver captures `beam_capture_fraction(link, radius_m, d)`,
`1 − exp(−2 r² / w²)`, of the beam through its aperture of radius `radius_m`.
It delivers that times its `efficiency` to its bus, up to its `max_power_w`
(`link_received_power_w`). With a 1.064 µm beam from a 5 cm waist, for
example, a receiver of 25 cm radius catches nearly all of the beam out to
about 10 km, about a quarter at 100 km, and under 1 % at 1000 km.

The delivered power is generation in the receiver's power budget. It shows in
`power_generated_w` and in the receiver's own `power_received_w` column,
charges its battery, and helps an empty battery recover. A full receiver
battery curtails it, as it does solar power, while the transmitter keeps paying
for the beam.

The link can close within `range_m`, while the segment between the spacecraft
clears the planet (a sphere of its equatorial radius). When the transmitter's
battery sheds tier 1, the beam is unavailable until the tier clears. A
`PowerReceiver` can only receive a power beam, and a `LaserTerminal` can only
transmit one (`terminal_supports(terminal, link, side)`).

Policies that score links by their forces never select a power beam. Use
`PowerDeliveryPolicy()`, which scores each link by the power it would deliver,
or `EnergyBalancingPolicy()`, which charges the lowest state of charge (see
[Energy balancing](#Energy-balancing)). To weigh momentum and power links
together, use a [weighted mix](#Weighted-mixes).

A run with a power beam adds the `power_received_w` and `power_beamed_w`
columns per spacecraft to the
[power columns](output_columns.md#Power-budget).

Pointing and jitter, the receiver's incidence angle, the atmosphere, beam
quality and the beam's radiation pressure are not modeled.

## Data links

A `DataLink` carries data one way, from a transmitter terminal to a receiver,
and applies no force. Its `budget` sets the rate. In this fragment `sc`, `bus`,
`other`, `model` and `spacecraft` are those of the examples above:

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

There are two built-in budgets, and you can write your own:

- `FixedRateBudget(; rate_bps)` gives a constant rate at any range up to
  `range_m`. It also works between virtual terminals.
- `ShannonBudget(; bandwidth_hz, frequency_hz, implementation_efficiency=0.5,
  max_rate_bps=Inf)` is a radio link budget between two `RadioTerminal`s. The
  received power follows Friis, `P_r = P_t G_t G_r (λ / 4πd)²`, with the
  link's `transmit_power_w` as `P_t`, the two radios' `gain_dbi`, and
  `λ = c / frequency_hz`. The noise is `k_B T bandwidth_hz`, with the
  receiver's `noise_temperature_k` as `T`. The rate is
  `implementation_efficiency × bandwidth_hz × log₂(1 + P_r / N)`, up to
  `max_rate_bps`. For the link above that is about 5.7 Mbit/s at 10 km,
  2.4 Mbit/s at 100 km and 0.17 Mbit/s at 1000 km.
- A budget kind of your own subtypes `AbstractDataBudget` and implements
  `budget_rate_bps(budget, link, transmitter, receiver, distance_m)`. The
  terminals are `nothing` at a virtual end.

The budget's rate is then capped at each `RadioTerminal`'s `max_rate_bps`
(`link_data_rate_bps`).

A `RadioTerminal` can transmit or receive. While it transmits, it runs at the
link's `transmit_power_w` and draws `transmit_power_w / efficiency` as a tier-1
load. A receiving radio draws nothing, so give its electronics a constant
`PowerLoad`. When the transmitter's battery sheds tier 1, the link is
unavailable until the tier clears. A `LaserTerminal` or a `PowerReceiver`
cannot carry a data link.

In a run with a data link, every spacecraft has two more integrated states,
`data_sent_bits` and `data_received_bits`. Their rates are the sums of the
rates its active links carry out and in, evaluated at every right-hand-side
call, so the solver integrates these counters with the other states under its
integration tolerances. They continue across a checkpoint resume.

A data link is available on the same terms as a power beam: within `range_m`,
while the segment between the spacecraft clears the planet. Two-way traffic is
two links, one each way. A radio carries one active link at a time, so of a
two-way pair between the same two radios the scheduler selects the direction
with the higher score.

Use `DataThroughputPolicy()`, which scores each link by the rate it would
carry; the force and power policies never select a data link. A run with a
data link adds the `data_rate_tx_bps`, `data_rate_rx_bps`, `data_sent_bits` and
`data_received_bits` columns per spacecraft (see
[the data columns](output_columns.md#Data-links)).

Pointing and antenna patterns, Doppler, atmospheric and polarization losses,
interference between links and acquisition time are not modeled. Onboard
storage and ground stations are covered next.

### Data storage

A spacecraft can hold data in a `DataStore` and produce it with data sources.
The data links then carry only what the stores hold, or what passes through
them. In this fragment `sc` and its root link `bus` are as above:

```julia
add_data!(sc, bus, DataStore(capacity_bits=8e9, initial_bits=0.0))   # the root link, at most one
add_data!(sc, bus, DataSource(rate_bps=2e5))                          # any link, any number
```

In a run with a store, every spacecraft has three more integrated states,
which stay zero for a spacecraft without a store: `data_stored_bits`, the data
its store holds; `data_produced_bits`, what its sources have produced; and
`data_lost_bits`, the production a full store could not keep. They balance,
`stored = initial + produced + received − sent − lost`, where
`data_sent_bits` and `data_received_bits` count what the links actually carry.

A link carries at most its rate (`link_data_rate_bps`), and less at a store,
so links lose nothing:

- An empty store passes on what comes in: its production and what its
  incoming links carry, shared among its outgoing links in proportion to their
  rates. A relay with an empty store therefore forwards data without holding
  any.
- A full store takes in no more than it sends on. Its production beyond that
  is lost and counted in `data_lost_bits`.
- A spacecraft without a store neither limits its links nor is limited by
  them.

The rates are found at every right-hand-side evaluation, along chains of
relays.

A store counts as full when it holds within a relative `1e-9` of its
capacity, and as empty when it holds no more than `1e-9` times its capacity.
Events mark it full or empty when it crosses either boundary, as a battery's
do, and the flags are checked again at the end of each step. The flags
continue across a checkpoint resume, and a resume with another capacity or
initial content is refused.

`DataSource(; rate_bps, name="")` produces a constant rate while its
spacecraft is active. A source kind of your own subtypes `AbstractDataSource`
and implements `data_production_bps(source, x::StateSample, ctx::LinkContext,
t)`, its rate in bit/s: for example imaging only over a region, a duty cycle,
or a rate that depends on the attitude. The engine may call it more than once
at the same state and time, so it must not advance any state of its own. A
run refuses at setup a source on a spacecraft without a store, and a source
kind without a method. As with a battery, a run with a data store needs a
solver mode that integrates the full right-hand side.

`DataBacklogPolicy(; horizon_s)` schedules by backlog. It scores a data link
by the data it could move in the next `horizon_s` seconds, as a rate in bit/s:

```text
min(c, stored_tx / horizon_s + production_tx, room_rx / horizon_s + outflow_rx)
```

Here `c` is the link's rate, `production_tx` the transmitter's current source
production rate, and `outflow_rx` the rate the receiver currently sends on.
The outflow is sampled once per scheduling step, before the new selection,
from the previously active links that are still available, with the same flow
limits the run applies; the policy does not assume the links it is about to
select will be active. Production is sampled at the current state and time
and assumed to continue over the horizon, so future changes of a source or of
incoming relay traffic are not predicted. The score is an estimate for
choosing links; the actual transfers still obey the stores' flow limits.

An empty sender that is still producing can keep its link active. A sender
with neither stored data nor production, or a receiver with neither free space
nor outflow, scores zero. A transmitter without a store counts as always
holding data, and a receiver without one as always having room. The policy
shares `DataThroughputPolicy`'s unit, so the two mix in a `WeightedPolicy`
without rescaling. A policy of your own reads the stores with
`stored_bits(ctx, i)` (NaN without a store) and `store_room_bits(ctx, i)` (Inf
without one).

In relay networks the backlog score has limits to keep in mind. It estimates a
local transfer, not additional data delivered to the final destination. A
full relay can attract incoming traffic that displaces its own production
while a better direct path goes unused, and two full stores can keep
exchanging data without increasing delivery anywhere. Each incoming candidate
receives the receiver's whole outflow credit, which is neither reduced by the
receiver's own production nor divided among the candidates. When several links
are selected, the actual transfers still share the available flow and respect
the store limits. Compare delivered and lost data when you choose a policy for
a relay network.

A run with a store adds the `data_stored_bits`, `data_store_fill` (the
fraction of the capacity held), `data_produced_bits` and `data_lost_bits`
columns per spacecraft (see
[the storage columns](output_columns.md#Data-storage)). `data_rate_tx_bps` and
`data_rate_rx_bps` then show the rates the links carry under the store
limits.

One store holds one stream of data. There are no priority classes, no queues
per destination, and no age or latency of the data held.

### Ground stations

A ground station is a fixed point on the rotating planet, and its antennas are
link terminals. In this fragment `sender` and `link` are the radio and the data
link of the example above:

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

`GroundStation(; name, lat_deg, lon_deg, alt_m=0.0, min_elevation_deg=10.0,
antennas)` needs a name of its own within the model, a latitude within
[-90, 90], a mask within [0, 90) degrees and at least one antenna.
`add_ground_station!` returns the station's number `s`, and station `s` is
the node `-s`: its antenna `k` is the endpoint `(-s, k)` in link keys, the
history, the contact windows and the link schedule. A link joins a spacecraft
and a station; two stations cannot link. `forbidden_pairs` may name `(i, -s)`.
`station_pairs` lists `(i, -s)` pairs, so the spacecraft transmits; for
uplinks, reverse each pair.

The station sits at its geodetic latitude, longitude and altitude on the
planet's ellipsoid, and turns with the planet. Its inertial position is
`L_PIᵀ r` and its velocity `L_PIᵀ (ω × r)`, with `r` its planet-fixed
position, `L_PI` the planet's frame and `ω` its spin. Availability and the
scheduling policies read the frame at each accepted step. The rates read it at
the time of each right-hand-side evaluation, so a station moves through a step
instead of holding its place at the step's start.

A link to a station is available while the spacecraft is at least
`min_elevation_deg` above the station's local horizon (the plane normal to its
geodetic up direction) and within the link's `range_m`. The mask stands for
the line of sight: a `DataLink` does not also test the segment to a station
against the planet.

The antennas are `RadioTerminal`s, numbered `1:length(antennas)`, and each
carries one active link at a time. They steer freely above the mask, so their
`boresight` is not read. A station has no power budget: its antennas draw
nothing and are never shed.

A downlink counts in its spacecraft's `data_sent_bits`, and an uplink, a link
from `(-s, k)` to a spacecraft's radio, in that spacecraft's
`data_received_bits`. In a run with data links each station also counts the
bits it has received and sent, in two integrated states
(`u.gs[s].received_bits` and `u.gs[s].sent_bits`) written as the
`gs{s}_received_bits` and `gs{s}_sent_bits` columns (see
[the data columns](output_columns.md#Data-links)). They continue across a
checkpoint resume.

With data stores, a station acts as a store with unbounded capacity and
supply, so it never limits a link: a downlink carries what its spacecraft's
store sends, and an uplink what its spacecraft's store can take in.
`stored_bits(ctx, -s)` and `store_room_bits(ctx, -s)` are Inf, so
`DataBacklogPolicy` drains the fullest stores to the ground first.

Station links are candidates like any other. A station is never an orbital
policy's target. It has no battery, so its `state_of_charge` is NaN and
`EnergyBalancingPolicy` scores a link to it 0.

In the link schedule a station end shows as `-s` in `a_spacecraft` or
`b_spacecraft`, and by its registered name in `a_name` or `b_name`; these name
columns appear only in runs with stations. Saved scenes include each station's
planet-fixed position. SpaceAGORAVisualizer 0.2 draws each station as a green
marker with its name on the rotating globe, together with its active and
contact links in either direction; the labels toggle shows or hides the names,
and the interlinks toggle shows or hides the stations and their links. A run
without stations writes no station entries in its scene and no name columns in
its schedule.

A checkpointed run with stations records each station's site, mask and
antennas, and refuses to resume when any of them has changed.

Terrain and horizon masks, antenna steering limits and patterns, and
atmospheric losses are not modeled.

## Scheduling

At the start of the run and after every accepted solver step, the engine marks
a candidate available when its spacecraft pair is not forbidden, both
spacecraft are active, neither terminal is shed, and its link's
`link_available` allows it. The run's scheduling policy then scores each
available link, and the scheduler selects the terminal matching with the
largest total score. Only links with a positive score can be selected. The
selection is held between accepted steps, so the solver's trial stages and
rejected steps never change it. An inactive endpoint contributes no force even
if an earlier selection still names its link.

In a run with interlinks, every spacecraft has two more integrated states.
`interlink_dv` integrates the acceleration that its selected links' forces
give (m/s, inertial). `interlink_delta_sma` integrates the rate at which those
forces change the osculating semimajor axis, `2a²(v·F)/(μm)` (m); it is not a
comparison against a separately propagated unforced trajectory. Both stay zero
for links that apply no force, such as power beams and data links. They are
saved in the [interlink columns](output_columns.md#Interlinks). The history
belongs to the configuration that ran, so read it from the returned run, as in
the examples above. With the default state isolation, the configuration you
passed in is left unchanged.

### Policies

A scheduling policy is a value of a policy type. Give the run one with
`SimulationConfiguration(scheduling_policy_model=policy)`; the default is
`SemimajorAxisRatePolicy()`. The run keeps it in a `SchedulingPolicyModel`,
the scheduler, which also holds the matching workspace, and
`args.scheduling_policy_model.objective` is the policy.

| Policy | Scores each available link by | Unit |
| --- | --- | --- |
| `SemimajorAxisRatePolicy(; targets=1)` | the rate at which its force raises the osculating semimajor axis of its target spacecraft, `2a²(v·F)/(μm)`, summed over the target ends | m/s |
| `EccentricityRatePolicy(; targets=1)` | the rate at which its force raises its target's eccentricity, `ê·ė` | 1/s |
| `PowerDeliveryPolicy()` | the electrical power it would deliver to the receiver's bus, `link_received_power_w` | W |
| `DataThroughputPolicy()` | the data rate it would carry, `link_data_rate_bps` | bit/s |
| `DataBacklogPolicy(; horizon_s)` | the data it could move in the next `horizon_s`: the least of its rate, the transmitter's stored data over `horizon_s` plus its current production rate, and the receiver's free space over `horizon_s` plus its current outflow | bit/s |
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
charge, `P × (s_tx − s_rx)`. As a result:

- power flows only from a fuller battery to an emptier one, and a beam into a
  full battery never scores. A battery the power system counts as full (its
  stored energy within a relative `1e-9` of its capacity) counts as 1, so two
  full batteries tie at 0;
- of receivers competing for one sender's terminal, the emptier one wins at
  equal delivered power, while a much better-coupled receiver can still win
  over a slightly emptier one;
- the score is signed: a beam from an emptier battery to a fuller one scores
  negative, so it is never selected alone and counts against the beam in a
  mix;
- a link with no battery at either end, or one that would deliver no power,
  scores 0.

The selection changes only at accepted steps, so when two states of charge
cross within a step, the beam runs to the step's end. `PowerDeliveryPolicy`,
by contrast, beams into a full battery: the receiver curtails, and the sender
still pays.

### Weighted mixes

`WeightedPolicy(PolicyTerm(policy; scale, weight=1.0)...)` scores each link by
`Σ weight × (score / scale)` over its terms, summed in term order. `scale` is
in the term's own unit, the amount of its score that counts as one, so each
term is dimensionless before the weights apply. It must be finite and
positive. `weight` must be finite and may be negative. The scales are fixed,
not normalized by each step's best link, so a link's score does not depend on
the other candidates.

```julia
policy = WeightedPolicy(
    PolicyTerm(SemimajorAxisRatePolicy(targets=1); scale=1.0),   # 1 m/s of SMA rate counts as one
    PolicyTerm(EnergyBalancingPolicy(); scale=100.0))            # and so do 100 W
args = SimulationConfiguration(...; interlink_model=model, scheduling_policy_model=policy)
```

Each term applies its own rules: an orbital term scores power beams 0 and a
power term scores momentum links 0, which is what lets a mix weigh link kinds
against each other. Each term's checks also run, and an eccentricity term
keeps its near-circular refusal over the mix's candidates. Mixes nest.

### A custom policy

A policy type subtypes `AbstractSchedulingPolicy` and implements
`candidate_score`. This one prefers short links; `...` and `model` are as in
the examples above:

```julia
struct ShortLinkPolicy <: SpaceAGORA.AbstractSchedulingPolicy end
SpaceAGORA.candidate_score(::ShortLinkPolicy, candidate, ctx, t) =
    1 / sqrt(sum(abs2, candidate.a.pos_ii - candidate.b.pos_ii))

args = SimulationConfiguration(...; interlink_model=model, scheduling_policy_model=ShortLinkPolicy())
```

The candidate (`LinkCandidate`) has its `key`, its `link`, its endpoints'
state samples `a` and `b` (in key order, so a directed link's transmitter is
`a`), and the terminal hardware at its ends, `terminal_a` and `terminal_b`
(`nothing` for a virtual terminal). The context (`SchedulingContext`) has the
planet's gravitational parameter `mu`, the `link_context` to pass to the link
hooks, and the run's `spacecraft`. `state_of_charge(ctx, i)` gives spacecraft
`i`'s battery state of charge (NaN without a battery), and `stored_bits(ctx,
i)` and `store_room_bits(ctx, i)` give the data its store holds and the room
left (NaN and Inf without a store). `t` is the elapsed time in seconds.

The hook is called for available candidates only; the others score 0. A
nonfinite score is refused, and only positive scores can be selected. A run
refuses at setup a policy type without a `candidate_score` method.

A policy holds no run state: it is a plain value whose fields say what it
scores. A checkpoint records its type and every field by value (numbers with
their bits; arrays, tuples, named tuples and structs element by element and
field by field), so an equal policy rebuilt for the resume resumes, and a
resume with another value is refused. A function field is recorded by its
type and, for a closure or a callable struct, every value it holds. A named
function is recorded by its type alone, so a change to its code is not
detected. A closure reads the globals it uses when it is called, and their
values are not recorded, so build the closure in a function or a `let` block
where it captures them. A checkpointed run refuses at setup what a checkpoint
cannot record by value:

- a closure that captures a variable reassigned after the closure was made;
- a dictionary or a set;
- a stream, a task or a lock;
- a pointer, or a value holding one (such as a `Regex`), whose address differs
  from one process to the next;
- an unassigned field or element;
- fields that contain themselves.

### Contact windows

The history (`model.history`) has one `InterLinkHistory` sample per accepted
step: its time, the links the scheduler selected (`active`) and the links that
were available (`available`). The available links are recorded in a run with a
data link, or when you ask for them with
`InterlinkModel(; record_available=true)`; `record_available=false` leaves
them out. In a large constellation over a long run they make a long list at
every step, which is why other runs leave them out by default.
`contact_windows(history)` turns the history into windows, one
`(key, start_s, end_s, duration_s)` for each span over which a link could
close; `contact_windows(history; active=true)` gives the spans over which it
was selected. `contact_time_s(history, key)` sums a link's windows.

```julia
solution = run_simulation(args; return_solution=true)
history = solution.prob.p.args.interlink_model.history
windows = contact_windows(history)
contact_time_s(history, ((1, 1), (2, 1)); active=true)
```

A window starts at the first accepted step that has the link and ends at the
next step that does not, or at the last step, so its resolution is the
scheduler's. A checkpointed run keeps its history in
`simulation_checkpoint.interlinks.bin` beside its checkpoint, and a resumed run
takes it back, so a resumed history holds the whole run, as the results do
(see [Outputs](outputs.md)). If that record is missing or unreadable, the run
warns and records only the resumed interval; a resume that solves nothing then
returns an empty schedule. With the same configuration, checkpoint schedule
and deterministic models whose complete run state is restored, a resumed
history matches the uninterrupted checkpointed run. A run without checkpoints
can take different adaptive steps and so record different scheduling samples.

### The link schedule

Each history sample also records the score that selected each active link
(`scores`), read from the scheduler after it selects. `link_schedule` turns the
history into a table with one row per window: a span over which a link was
selected (`window == "active"`) and, when the history records the available
links, a span over which a link could close (`"contact"`). Each row gives the
link (its number, kind, direction and endpoints), the span (`start_s`,
`end_s`, `duration_s` and the number of samples, `steps`) and, for an active
row, the score at its first and last sample and its extremes.

```julia
results = run_simulation(args; return_results=true)
schedule = link_schedule(results)              # or link_schedule(model; policy=scheduler)
active = schedule[schedule.window .== "active", :]
```

The scores are in the policy's unit (see [Policies](#Policies)): m/s for
`SemimajorAxisRatePolicy`, 1/s for `EccentricityRatePolicy`, W for the power
policies, bit/s for the data policies and 1 for a `WeightedPolicy`. The
table's metadata names the policy and the unit. With
`InterlinkModel(; record_terms=true)` and a `WeightedPolicy`, each sample also
records each top-level term's own score, in its own unit, and the table gains
`term1_start`, `term2_start`, ...: each term's score at the span's first
sample. A run with another policy refuses `record_terms=true`. Recording the
terms leaves the mix's scores, and so the schedule, unchanged.

A run that writes results writes the table beside them, as
`simulation_results_link_schedule.csv` (when `save_csv` is set) and
`simulation_results_link_schedule.feather`, and lists both in the results
manifest with the policy and its unit (see the
[link schedule's columns](output_columns.md#Link-schedule)). This also holds
for long-format results and results streamed at each checkpoint. The table is
built from the history, so thinning the results rows does not thin it. A run
without interlinks writes none, and `link_schedule(results)` returns the
table's columns with no rows. Like the contact windows, the schedule has the
scheduler's resolution. With `save_visualization_scene`, the viewer payload
carries the schedule's windows, so the viewer can draw each link while it is
active (see [the viewer page](outputs.md#Interactive-viewer-page)).

A resumed run restores its recorded history from beside the checkpoint,
including earlier windows, and its metadata and manifest give
`history_start_s` and `history_end_s`. If the history is unavailable, a warning
names the missing earlier interval; read that interval as not recorded, not as
an interval without links. A few resume cases need care:

- A resume that only finishes interrupted trajectory file writes rebuilds its
  schedule from the validated checkpoint history. Without that history, it
  keeps any earlier schedule file but does not list it in the recovered
  manifest or the returned files.
- A completed resume that solves nothing returns the validated history
  without rewriting files.
- An unusable history file is kept under a `before_resume_at_<t>s` name, just
  before the first replacement checkpoint is written. A resume that fails
  before then leaves it untouched.

### Optional matching backend

Tree-shaped components of the candidate graph are solved in Julia. Of parallel
links (directed links both ways between the same two terminals), only the
higher-scoring one is a candidate, so a two-way pair alone is a tree. Cyclic
components use NetworkX's maximum-weight matching through the optional
PythonCall extension. Loading the core package and tree matching do not start
Python. A cyclic component without the backend raises an error that says what
to install, rather than falling back to an approximate match.

To install the backend, add PythonCall to the Julia environment used for the
run, and install NetworkX (`>=3.4,<4`) in a compatible Python environment. Set
`JULIA_PYTHONCALL_EXE` to that interpreter before starting Julia; when you
manage Python yourself, also set `JULIA_CONDAPKG_BACKEND=Null`. Then load
`PythonCall` and `SpaceAGORA` explicitly, in either order. For example, in a
separate Julia project containing SpaceAGORA:

```julia
using Pkg
Pkg.add("PythonCall")
using SpaceAGORA, PythonCall
```

The backend needs Python's global interpreter lock. Run cyclic matching on the
thread that owns it, or use independent process workers; adding Julia threads
does not make arbitrary calls into Python safe. Both matching paths use exact
integer weights derived from the floating-point scores. When scores tie, the
two paths, or two backend versions, may choose different matchings of equal
total score.

## A custom link kind

A kind implements `link_range_m` and the effects it has. This one pulls its
endpoints together with a constant force; `model` and `spacecraft` are as in
the examples above:

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

A run refuses at setup a link whose kind has no `link_range_m` method. Links
of different kinds can be registered in one model.

At a ground station end (see [Ground stations](#Ground-stations)), the hooks
receive the station's state sample: its `spacecraft` is the `GroundStation`,
its `spacecraft_index` is `-s`, and it has no mass and no attitude. The station
sits on the planet, so a kind that tests occultation by the planet should
leave a station end to the station's elevation mask, as `DataLink` does
(`a.spacecraft isa GroundStation`). A station's antennas are `RadioTerminal`s,
which carry a kind only where `terminal_supports` allows it.

### Terminals and torques

A kind that reads its terminal hardware, or that turns the spacecraft it
links, implements the terminal-aware forms of the hooks. The engine calls these
forms, and their defaults call the four-argument hooks above, so a kind that
defines only those gives the same result, bit for bit.

- `link_available(link, terminal_a, terminal_b, a, b, ctx)` decides whether the
  link can close, given the terminal hardware at its two ends (`nothing` for a
  virtual terminal).
- `interlink_wrench(link, terminal_a, terminal_b, a, b, ctx)` returns
  `(force_a_ii, torque_a_body, force_b_ii, torque_b_body)`. The forces act as
  `link_force`'s do, and the scheduling policies score candidate links by
  them. Each torque is in its spacecraft's body frame, about its body reference
  point, as a mounted actuator's torque is (see
  [Frames](../reference/extension_conventions.md#Frames)). The engine adds each
  torque to its spacecraft's angular acceleration as `I⁻¹ τ`, with `I` the
  inertia tensor, which is exact for the rigid spacecraft that interlinks
  support.
- `link_applies_torque(link) = true` says that a kind's torques can be
  nonzero. A run with such a link needs attitude (`orientation_sim=true`), and
  setup refuses it otherwise. A nonzero torque from a kind that does not
  declare this stops a run without attitude.
- `terminal_pose(spacecraft, terminal)` gives a mounted terminal's boresight
  and location in its spacecraft's body frame. A force applied at that
  location has the torque the terminal's own recoil would have: the cross
  product of the location and the force, in the body frame.

### Supported scope

Interlinks apply forces, and through `interlink_wrench` torques, to rigid
spacecraft, or carry power or data. The laser momentum link does not model
terminal pointing, occultation, beam losses or torques; `LaserBeamLink` adds
them as options (see [Laser physics](#Laser-physics)). Neither models the
terminals' thermal dynamics. Links through hardware terminals draw battery
power and close when it is shed; links between virtual terminals draw none.
Actuator faults on a link-driven terminal do not act on its link. The
spacecraft's `battery_energy_index` and `tempurature_index` fields do not
constrain link availability (`battery_energy_index` only sets a battery's
default initial state of charge). The cavity magnification and the force law
are modeling assumptions; passing the implementation tests does not validate
a physical laser system.

Use a full-state solver such as `:tsit5`, or the opt-in `:multirate_mri` route
described in
[Solver configuration](solver_configuration.md#Coupled-multirate-integration).
A run with interlinks refuses the split modes `:split_imex`, `:multirate` and
`:gravity_backbone_split`. It also refuses linked spacecraft that are
externally propagated, articulated or compliant, and
`run_constellation_ensemble` refuses an interlinked constellation. Near
circularity (an eccentricity of at most `sqrt(eps(Float64))`),
`EccentricityRatePolicy` also refuses candidate sets that could apply several
positive-score links to one target at once, because the eccentricity-rate
norms of individual links do not add there. `SemimajorAxisRatePolicy` has no
such restriction.

Scheduling follows accepted steps, so changes of range or of the preferred
link are resolved at that cadence. Check timestep convergence for a new case.

A checkpoint resume must keep the interlink topology, the link kinds (each
named by its module and type, with its concrete type parameters) and their
parameters, the scheduling policy (its type and every field), and each
spacecraft's terminal hardware: its kinds, in mounting order, and every field
of each, such as a receiver's aperture or a transmitter's efficiency. Start a
new run when any of these change. Checkpoints written by SpaceAGORA 0.4 and
earlier still resume, and so do checkpoints written with the symbol policies
of 0.5, under the typed policy that schedules the same. A resumed run takes its
interlink history back from beside its checkpoint, as its results take back
their rows, so its history, schedule and results table hold the whole run,
and the integrated diagnostic states continue from the checkpoint. Restoring
the recorded history requires the same `record_available` and `record_terms`
settings; a resume that changes either is refused before any result stream is
changed. Start a new run to change these settings.

## Deprecated names

The names from when every link was a laser link keep working until SpaceAGORA
1.0. The old keywords, fields and functions warn once per process, and throw
under `julia --depwarn=error`. The two old type names are aliases of the new
types and cannot warn on their own: `InterLinkParameters()`,
`InterLinkParameters(terminal)`, the positional
`InterLinkParameters(P, B, range)` and
`ConstellationInterLinkModel(; forbidden_pairs)` run silently, so search your
code for those names when you move to the new ones. The `laser_dv` and
`laser_delta_sma` columns cannot warn either: a run writes them beside the
`interlink_dv` and `interlink_delta_sma` columns, with the same values, until
1.0. Leave the `:laser_dv` and `:laser_delta_sma` fields out of an explicit
`save_fields` list to stop writing them.

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
without a deprecation period, together with their `target_index` and
`target_indices` keywords and properties and the older `target_idx` and
`target_idxs`. A call with a symbol throws an `ArgumentError` that names its
typed replacement, with or without `julia --depwarn=error`, and so does
reading or setting a removed property, which names `.objective` or
`.objective.targets`. Checkpoints written with the symbol policies still
resume under the typed policy that schedules the same.

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

The two integrated states the interlinks add to each spacecraft were renamed
in 0.6.0 as well. The name of a state component cannot warn, so they were
renamed without a deprecation period; their columns keep the old names as
aliases until 1.0 (see [Deprecated names](#Deprecated-names)). Checkpoints keep
the old names inside, so checkpoints written before 0.6.0 resume.

| Removed | Use instead |
| --- | --- |
| the state `laser_dv`, as in `solution.u[end].sc[i].laser_dv` | `interlink_dv` |
| the state `laser_delta_sma` | `interlink_delta_sma` |
