# The link contract: what one registered interlink does between its two terminals.
# Each registered link is one value of a link kind; the engine asks it whether it
# can close and what it does to its two endpoints. Included by `InterLinkModels`
# (interlink_graph.jl).

"""
    AbstractInterlink

Supertype of the link kinds. One value describes one registered link between two
terminals: register it with `register_link!`. A kind implements `link_range_m`,
and the effects it has: `link_force` for a link that pushes its endpoints (the
default exerts none), or `interlink_wrench` for one that also turns them.
`link_available` decides when the link can close (by default, within
`link_range_m`). The built-in kinds are `LaserMomentumLink`, `LaserBeamLink`,
`PowerBeamLink` and `DataLink`.
"""
abstract type AbstractInterlink end

"""
    LinkContext(; time_s=NaN, planet=nothing)

What a link hook may read besides its two endpoints: the time (s into the run)
and the run's planet. The engine passes the run's values; a call outside a run
gets the defaults.
"""
struct LinkContext{P}
    time_s::Float64
    planet::P
end
LinkContext(; time_s::Real=NaN, planet=nothing) = LinkContext(Float64(time_s), planet)

"""
    link_range_m(link) -> Float64

The maximum separation in metres at which `link` can close. Every link kind
implements it.
"""
function link_range_m end

"""
    link_available(link, a::StateSample, b::StateSample, ctx::LinkContext) -> Bool
    link_available(link, terminal_a, terminal_b, a::StateSample, b::StateSample, ctx::LinkContext) -> Bool

Whether `link` can close between its endpoints `a` and `b`, in the order of the
link's key (the lower `(spacecraft, terminal)` pair first). The engine asks at the start of the run and after every accepted step,
once the run's own rules (forbidden pairs, active spacecraft) allow the pair.

- The engine calls the six-argument form, with `terminal_a` and `terminal_b`
  the terminal hardware at the two ends (`nothing` for a virtual terminal). A
  kind that reads its terminals, such as their pointing (`terminal_pose`),
  implements it.
- Default: the six-argument form calls the four-argument form, which is true
  while the endpoints are apart and within `link_range_m(link)`.
"""
function link_available(link, a, b, ctx)
    separation = norm(a.pos_ii - b.pos_ii)
    return 0 < separation <= link_range_m(link)
end
link_available(link, terminal_a, terminal_b, a, b, ctx) = link_available(link, a, b, ctx)

"""
    link_directed(link) -> Bool

Whether `link` runs one way, from its first endpoint (the transmitter) to its
second (the receiver). A directed link is keyed in the order its endpoints were
registered, so the hooks receive the transmitter as `a`. An undirected link is
keyed with the lower `(spacecraft, terminal)` first. Default: undirected, as
for `LaserMomentumLink`.
"""
link_directed(link) = false

"""
    link_terminal_levels(link) -> (level_a, level_b)

The levels at which the hardware terminals at `link`'s two endpoints run while
it is active, in key order and in their actuators' units (watts of transmitted
optical power for a `LaserTerminal`). Each such terminal draws
`actuator_consumption(terminal, level, ctx).power_w` from its spacecraft's power
budget, and is at level 0 while its link is inactive. A virtual terminal draws
nothing. Default: zero levels.
"""
link_terminal_levels(link) = (0.0, 0.0)

"""
    link_terminal(actuator) -> Bool

Whether a mounted actuator is interlink terminal hardware, which can be a link
endpoint. True for `LaserTerminal`, `PowerReceiver` and `RadioTerminal`; a terminal kind of your own opts in with
`SpaceAGORA.link_terminal(::MyTerminal) = true`. A spacecraft's terminals are
its `n_terminal` virtual terminals, numbered `1:n_terminal`, followed by its
mounted terminals in mounting order. Default: false.
"""
link_terminal(actuator) = false
link_terminal(::LaserTerminal) = true

"""
    terminal_supports(terminal, link) -> Bool
    terminal_supports(terminal, link, side) -> Bool

Whether the terminal hardware `terminal` (a mounted actuator with
`link_terminal`) can carry `link`. `register_link!` refuses a link that either
endpoint's terminal cannot carry. A virtual terminal is refused when the link
requires hardware at both ends: a `DataLink` with a `ShannonBudget`, or a
`LaserBeamLink` with mounted recoil or terminal pointing.

- `side` is the end the terminal would take: `:transmitter` or `:receiver` for a
  directed link, `:either` for an undirected one. A kind whose terminals work
  one way implements the three-argument form.
- Default: the two-argument form, which is true. A `LaserTerminal` carries
  `LaserMomentumLink` and transmits a `PowerBeamLink`; a `PowerReceiver`
  receives a `PowerBeamLink`; a `RadioTerminal` transmits or receives a
  `DataLink`. Other kinds need a method of their own.
"""
terminal_supports(terminal, link) = true
terminal_supports(terminal, link, side::Symbol) = terminal_supports(terminal, link)
terminal_supports(::LaserTerminal, link) = false

"""
    terminal_shed_tier(actuator) -> Int

The battery shedding tier that closes the links of the terminal hardware
`actuator` (a mounted actuator with `link_terminal`) while a link drives it: 1 sheds it with the laser terminals, 2 with
the torque rods, 3 with the reaction wheels, 0 never. A shed terminal draws
nothing, and its links are unavailable until the tier clears. Its demand counts
in that tier of the power budget (0: with the kinds never shed). 1 for
`LaserTerminal` and `RadioTerminal`. Default: 0.
"""
terminal_shed_tier(actuator) = 0
terminal_shed_tier(::LaserTerminal) = 1

"""
    terminal_pose(spacecraft, terminal) -> (boresight_body, location_body)

Where `terminal`, terminal hardware mounted on `spacecraft` (a
`SpacecraftModel`), points and sits, in the spacecraft's body frame: its
`boresight` (a `PowerReceiver`'s `normal`) and its `location`, carried from its
link's frame by the link's configured pose.

- `location_body` is the lever arm about the body reference point that the
  terminal's own recoil uses (`actuator_wrench`), so a link that applies its
  force there (`interlink_wrench`) has the terminal's torque.
- A terminal without a `location` sits at its link's origin. One without a
  `boresight` or `normal` has no pointing axis, and is refused.
- Two mounted terminals with equal fields are the same value. If they sit at
  different poses, they cannot be told apart and are refused; give them
  distinct fields.
"""
function terminal_pose end

# Whether `link` reads the terminal hardware at both its ends, so that a virtual
# terminal cannot carry it (a `DataLink` with a `ShannonBudget`), and why, for
# the refusal.
_link_needs_terminals(link) = false
_needs_terminals_reason(link) = ""

# Whether `link` reads its endpoints' attitude, so that a run without attitude
# (`orientation_sim`) refuses it at setup.
_link_needs_attitude(link) = false

"""
    terminal_count(spacecraft) -> Int

The number of interlink terminals on `spacecraft`: its `n_terminal` virtual
terminals and, on a `SpacecraftModel`, its mounted terminal hardware (see
`link_terminal`).
"""
terminal_count(spacecraft) = Int(spacecraft.n_terminal)

# The hardware terminal number `number` stands for on `spacecraft`: `nothing` for
# a virtual terminal. `SpacecraftModel` methods are in InterlinkTerminals.
_terminal_device(spacecraft, number::Integer) = nothing
# Every terminal's hardware on `spacecraft`, by terminal number, for repeated
# lookups: entry `number` is `_terminal_device(spacecraft, number)`, and
# `nothing` stands for a spacecraft without hardware terminals.
_terminal_devices(spacecraft) = nothing

# The capacity (Wh) of `spacecraft`'s battery, NaN without one: the denominator of
# a scheduling policy's state of charge. `SpacecraftModel` methods are in
# InterlinkTerminals, which comes after the power models.
_battery_capacity_wh(spacecraft) = NaN
# Whether `spacecraft`'s battery counts as full holding `energy_wh` (Wh), as the
# power system decides it; false without one.
_battery_within_full(spacecraft, energy_wh::Real) = false
# The capacity (bit) of `spacecraft`'s data store, NaN without one: the room a
# scheduling policy reads. `SpacecraftModel` methods are in InterlinkTerminals,
# which comes after the data parts.
_data_store_capacity_bits(spacecraft) = NaN
# Current source production at this endpoint. The concrete spacecraft method is
# supplied after data parts load, alongside the store-capacity methods.
_data_store_production_bps(spacecraft, sample, context, t::Float64) = 0.0

# Implemented with store flow control after the spacecraft model is loaded.
function _data_store_outflows_bps end
# Whether `spacecraft`'s data store counts as empty, or as full, holding `bits`,
# as the store's own events decide it; false without a store.
_data_store_within_empty(spacecraft, bits::Real) = false
_data_store_within_full(spacecraft, bits::Real) = false

# The terminal number of an endpoint's terminal: a number as given, or a mounted
# terminal object (resolved on a `SpacecraftModel` in InterlinkTerminals).
_terminal_number(spacecraft, terminal::Integer) = Int(terminal)
_terminal_number(spacecraft, terminal) = throw(ArgumentError(
    "A link endpoint's terminal is a terminal number, or a terminal mounted on a SpacecraftModel; got $(typeof(terminal))."))

const _ZERO_FORCE = SVector{3, Float64}(0.0, 0.0, 0.0)

"""
    link_force(link, a::StateSample, b::StateSample, ctx::LinkContext) -> (force_a_ii, force_b_ii)

The inertial forces (N) an active `link` puts on its endpoints `a` and `b`, in
the order of the link's key. The engine reads them through `interlink_wrench`, whose
default returns them with zero torques: it adds them to the endpoints'
accelerations at every right-hand-side evaluation, and the scheduling policies
score candidate links by them. Default: no force.
"""
link_force(link, a, b, ctx) = (_ZERO_FORCE, _ZERO_FORCE)

"""
    interlink_wrench(link, terminal_a, terminal_b, a::StateSample, b::StateSample, ctx::LinkContext) -> (force_a_ii, torque_a_body, force_b_ii, torque_b_body)

What an active `link` does to its endpoints `a` and `b`, in the order of the
link's key: the inertial force (N) on each, and the torque (N·m) on each in its
body frame, about its body reference point, as a mounted actuator's torque.
`terminal_a` and `terminal_b` are the terminal hardware at the two ends
(`nothing` for a virtual terminal).

- The engine adds the forces to the endpoints' accelerations at every
  right-hand-side evaluation, and the scheduling policies score candidate links
  by them. It adds each nonzero torque to its spacecraft's angular acceleration
  as `I⁻¹ τ`, with `I` the spacecraft's inertia tensor.
- A kind whose torques can be nonzero says so with `link_applies_torque`, so
  that a run without attitude (`orientation_sim`) refuses it at setup.
- Default: `link_force`'s two forces, with zero torques.
"""
function interlink_wrench(link, terminal_a, terminal_b, a, b, ctx)
    force_a, force_b = link_force(link, a, b, ctx)
    return (force_a, _ZERO_FORCE, force_b, _ZERO_FORCE)
end

"""
    link_applies_torque(link) -> Bool

Whether `interlink_wrench` can return nonzero torques for `link`. A run with such a
link needs attitude (`orientation_sim`): setup refuses it otherwise, and a
nonzero torque from a kind that does not say so stops a run without attitude.
Default: false.
"""
link_applies_torque(link) = false

# The state sample of a link endpoint: spacecraft `index` of a run state
# `view`, with its model when known. A run state holds every field; a velocity
# or mass a caller's partial state lacks reads NaN.
const _NAN3 = SVector{3, Float64}(NaN, NaN, NaN)
function link_endpoint(view, spacecraft, index::Integer)
    vel = hasproperty(view, :vel) ? SVector{3, Float64}(view.vel) : _NAN3
    mass = hasproperty(view, :mass) ? Float64(view.mass) : NaN
    q = hasproperty(view, :q) ? SVector{4, Float64}(view.q) : nothing
    ω = hasproperty(view, :ω) ? SVector{3, Float64}(view.ω) : nothing
    return StateSample(SVector{3, Float64}(view.pos), vel, mass, q, ω, spacecraft, Int(index))
end

# ---- The laser momentum link ---------------------------------------------------

# The speed of light, the one the laser terminal's recoil uses.
const SPEED_OF_LIGHT_MPS = _SPEED_OF_LIGHT_MPS

"""
    LaserMomentumLink(; power_w=10_000.0, magnification=100.0, range_m=200e3)
    LaserMomentumLink(terminal::LaserTerminal)

An open-cavity laser link that pushes its two spacecraft apart: an idealized
central force of `magnification * power_w / c` on each endpoint, directed away
from the other, equal and opposite, at any separation up to `range_m`. It models
no pointing, occultation or beam losses; `LaserBeamLink` adds them as options.
`power_w` (W) must be nonnegative, `magnification` and `range_m` (m) positive.
The terminal form takes a link
between two terminals of that design: its `max_power_w`, `magnification` and
`range_m`.
"""
struct LaserMomentumLink <: AbstractInterlink
    power_w::Float64
    magnification::Float64
    range_m::Float64

    function LaserMomentumLink(power_w::Real, magnification::Real, range_m::Real)
        isfinite(power_w) && power_w >= 0 || throw(ArgumentError("power_w must be finite and nonnegative."))
        isfinite(magnification) && magnification > 0 || throw(ArgumentError("magnification must be finite and positive."))
        isfinite(range_m) && range_m > 0 || throw(ArgumentError("range_m must be finite and positive."))
        new(Float64(power_w), Float64(magnification), Float64(range_m))
    end
end

function LaserMomentumLink(; power_w=_NOT_PASSED, magnification=_NOT_PASSED, range_m=_NOT_PASSED,
                           P=_NOT_PASSED, B=_NOT_PASSED, range=_NOT_PASSED)
    power = _laser_keyword(:P, P, :power_w, power_w, 10_000.0)
    gain = _laser_keyword(:B, B, :magnification, magnification, 100.0)
    reach = _laser_keyword(:range, range, :range_m, range_m, 200e3)
    return LaserMomentumLink(power, gain, reach)
end

LaserMomentumLink(terminal::LaserTerminal) =
    LaserMomentumLink(terminal.max_power_w, terminal.magnification, terminal.range_m)

link_range_m(link::LaserMomentumLink) = link.range_m

function link_force(link::LaserMomentumLink, a, b, ctx)
    separation = a.pos_ii - b.pos_ii
    distance = norm(separation)
    distance == 0 && return (zero(separation), zero(separation))
    force = (link.magnification * link.power_w / SPEED_OF_LIGHT_MPS) * (separation / distance)
    return (force, -force)
end

# Between laser terminals, both run at the link's power while it is active.
link_terminal_levels(link::LaserMomentumLink) = (link.power_w, link.power_w)
terminal_supports(::LaserTerminal, ::LaserMomentumLink) = true

# The `InterLinkParameters` names of the fields, P, B and range: deprecated.
const _LASER_OLD_FIELDS = (P=:power_w, B=:magnification, range=:range_m)

function _laser_keyword(old::Symbol, old_value, new::Symbol, new_value, default)
    old_value === _NOT_PASSED && return new_value === _NOT_PASSED ? default : new_value
    new_value === _NOT_PASSED || throw(ArgumentError(
        "LaserMomentumLink: pass $(new) only; $(old) is its deprecated name."))
    _interlink_deprecation((LaserMomentumLink, old),
        "The InterLinkParameters keyword $(old) is deprecated and will be removed in SpaceAGORA 1.0; " *
        "use LaserMomentumLink($(new)=...) instead: only the name changes.")
    return old_value
end

@inline function Base.getproperty(link::LaserMomentumLink, name::Symbol)
    if name === :P || name === :B || name === :range
        new = _LASER_OLD_FIELDS[name]
        _interlink_deprecation((LaserMomentumLink, name),
            "LaserMomentumLink.$(name) (InterLinkParameters.$(name)) is deprecated and will be removed in " *
            "SpaceAGORA 1.0; use .$(new) instead: only the name changes.")
        return getfield(link, new)
    end
    return getfield(link, name)
end

Base.propertynames(::LaserMomentumLink, private::Bool=false) = (:power_w, :magnification, :range_m)

# Deprecated interlink names warn once per process and name, and throw under
# `julia --depwarn=error`.
@noinline _interlink_deprecation(key::Tuple, msg::String) = _hook_deprecation(key, msg)
