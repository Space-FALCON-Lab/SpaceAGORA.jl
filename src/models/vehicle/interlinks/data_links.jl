# Data links: a one-way link that carries data at a rate its link budget sets,
# from a transmitter terminal (a `RadioTerminal`) to a receiver; and the contact
# windows of a run's links. Included by `InterLinkModels` (interlink_graph.jl),
# after links.jl and power_links.jl.

"""
    link_data_rate_bps(link, transmitter, receiver, a::StateSample, b::StateSample, ctx::LinkContext) -> Float64

The data rate (bits per second) that `link` carries from its transmitting end
`a` to its receiving end `b` while it is active. `transmitter` and `receiver`
are the terminal hardware at the two ends, `nothing` for a virtual terminal.
The engine integrates it, at every right-hand-side evaluation, into the sending
spacecraft's `data_sent_bits` and the receiving one's `data_received_bits`, and
`DataThroughputPolicy` and `DataBacklogPolicy` score candidate links by it. Default: 0.
"""
link_data_rate_bps(link, transmitter, receiver, a, b, ctx) = 0.0

"""
    AbstractDataBudget

Supertype of the link budgets that set a `DataLink`'s rate. A budget kind is an
immutable struct that implements `budget_rate_bps`. The built-in kinds are
`FixedRateBudget` and `ShannonBudget`.
"""
abstract type AbstractDataBudget end

"""
    budget_rate_bps(budget, link, transmitter, receiver, distance_m) -> Float64

The data rate (bits per second) that `budget` gives `link`, a `DataLink`, at a
separation of `distance_m` between its ends. `transmitter` and `receiver` are
the terminal hardware at the two ends, `nothing` for a virtual terminal. The
`DataLink` then caps it at each `RadioTerminal`'s `max_rate_bps`. Every budget
kind implements it.
"""
function budget_rate_bps end

"""
    FixedRateBudget(; rate_bps)

A link budget with a constant rate: `rate_bps` bits per second while the link is
active, at any range up to the link's `range_m`. It works between virtual
terminals as well as between radios.
"""
struct FixedRateBudget <: AbstractDataBudget
    rate_bps::Float64

    function FixedRateBudget(rate_bps::Real)
        isfinite(rate_bps) && rate_bps > 0 || throw(ArgumentError("rate_bps must be finite and positive."))
        new(Float64(rate_bps))
    end
end
FixedRateBudget(; rate_bps::Real) = FixedRateBudget(rate_bps)

budget_rate_bps(budget::FixedRateBudget, link, transmitter, receiver, distance_m) = budget.rate_bps

"""
    ShannonBudget(; bandwidth_hz, frequency_hz, implementation_efficiency=0.5, max_rate_bps=Inf)

A radio link budget: the Shannon capacity of the received signal over thermal
noise, scaled by `implementation_efficiency` (the share of capacity the modem
and coding reach, in (0, 1]) and capped at `max_rate_bps`. With the link's
`transmit_power_w` as `P_t`, both radios' gains `G_t` and `G_r` (from their
`gain_dbi`) and the receiver's `noise_temperature_k` as `T`:

- received power (Friis): `P_r = P_t G_t G_r (λ / 4πd)²`, with `λ = c / frequency_hz`;
- noise: `N = k_B T bandwidth_hz`;
- rate: `implementation_efficiency × bandwidth_hz × log₂(1 + P_r / N)`.

Both ends must be `RadioTerminal`s; `register_link!` refuses a virtual
terminal. It does not model pointing, antenna patterns, atmospheric or
polarization losses, Doppler or interference.
"""
struct ShannonBudget <: AbstractDataBudget
    bandwidth_hz::Float64
    frequency_hz::Float64
    implementation_efficiency::Float64
    max_rate_bps::Float64

    function ShannonBudget(bandwidth_hz::Real, frequency_hz::Real, implementation_efficiency::Real, max_rate_bps::Real)
        isfinite(bandwidth_hz) && bandwidth_hz > 0 || throw(ArgumentError("bandwidth_hz must be finite and positive."))
        isfinite(frequency_hz) && frequency_hz > 0 || throw(ArgumentError("frequency_hz must be finite and positive."))
        isfinite(implementation_efficiency) && 0 < implementation_efficiency <= 1 ||
            throw(ArgumentError("implementation_efficiency must be in (0, 1]."))
        max_rate_bps > 0 && !isnan(max_rate_bps) || throw(ArgumentError("max_rate_bps must be positive."))
        new(Float64(bandwidth_hz), Float64(frequency_hz), Float64(implementation_efficiency), Float64(max_rate_bps))
    end
end
ShannonBudget(; bandwidth_hz::Real, frequency_hz::Real, implementation_efficiency::Real=0.5, max_rate_bps::Real=Inf) =
    ShannonBudget(bandwidth_hz, frequency_hz, implementation_efficiency, max_rate_bps)

# Boltzmann's constant, J/K.
const BOLTZMANN_J_K = 1.380649e-23

function budget_rate_bps(budget::ShannonBudget, link, transmitter::RadioTerminal, receiver::RadioTerminal, distance_m)
    wavelength = SPEED_OF_LIGHT_MPS / budget.frequency_hz
    gains = exp10((transmitter.gain_dbi + receiver.gain_dbi) / 10)
    received_w = link.transmit_power_w * gains * (wavelength / (4pi * Float64(distance_m)))^2
    noise_w = BOLTZMANN_J_K * receiver.noise_temperature_k * budget.bandwidth_hz
    rate = budget.implementation_efficiency * budget.bandwidth_hz * log2(1.0 + received_w / noise_w)
    return min(rate, budget.max_rate_bps)
end
budget_rate_bps(::ShannonBudget, link, transmitter, receiver, distance_m) = throw(ArgumentError(
    "A ShannonBudget needs a RadioTerminal at both ends of its DataLink, got $(nameof(typeof(transmitter))) and " *
    "$(nameof(typeof(receiver)))."))

# Whether a budget reads the terminal hardware at both ends, so a virtual
# terminal cannot carry its link.
_budget_needs_terminals(budget) = false
_budget_needs_terminals(::ShannonBudget) = true

"""
    DataLink(; transmit_power_w, budget, range_m)

A one-way data link from a transmitter terminal to a receiver: a directed link
that carries data and applies no force.

- **The rate.** `budget` (a `FixedRateBudget`, a `ShannonBudget`, or a budget
  kind of your own) sets it at the current separation. It is capped at each
  end's `max_rate_bps` when that end is a `RadioTerminal`. While the link is
  active, the transmitter's spacecraft integrates it into `data_sent_bits`, and
  the receiver's into `data_received_bits`.
- **The transmitter.** A `RadioTerminal` transmitter runs at `transmit_power_w`
  radio watts while the link is active and draws `transmit_power_w / efficiency`
  as a tier-1 load: when its battery sheds tier 1, the link closes. A virtual
  terminal draws nothing.
- **Availability.** The link can close within `range_m`, while the segment
  between the two spacecraft clears the planet (a sphere of its equatorial
  radius). At a ground station end, the station's elevation mask stands for
  that test.
- **Direction.** Two-way traffic is two links, one each way. A radio carries one
  active link at a time.
"""
struct DataLink{B <: AbstractDataBudget} <: AbstractInterlink
    transmit_power_w::Float64
    budget::B
    range_m::Float64

    function DataLink(transmit_power_w::Real, budget::B, range_m::Real) where {B <: AbstractDataBudget}
        isfinite(transmit_power_w) && transmit_power_w >= 0 ||
            throw(ArgumentError("transmit_power_w must be finite and nonnegative."))
        isfinite(range_m) && range_m > 0 || throw(ArgumentError("range_m must be finite and positive."))
        new{B}(Float64(transmit_power_w), budget, Float64(range_m))
    end
end

DataLink(; transmit_power_w::Real, budget::AbstractDataBudget, range_m::Real) = DataLink(transmit_power_w, budget, range_m)

link_range_m(link::DataLink) = link.range_m
link_directed(::DataLink) = true
link_terminal_levels(link::DataLink) = (link.transmit_power_w, 0.0)

terminal_supports(::RadioTerminal, link) = false
terminal_supports(::RadioTerminal, ::DataLink) = true
link_terminal(::RadioTerminal) = true
terminal_shed_tier(::RadioTerminal) = 1

_link_needs_terminals(link::DataLink) = _budget_needs_terminals(link.budget)
_needs_terminals_reason(::DataLink) = " (a ShannonBudget reads both radios)"

# A ground station end sits on the planet, inside the occulting sphere: its
# elevation mask, which the graph checks, is the link's line of sight.
function link_available(link::DataLink, a, b, ctx)
    separation = norm(a.pos_ii - b.pos_ii)
    0 < separation <= link.range_m || return false
    (_station_end(a) || _station_end(b)) && return true
    return !_segment_occulted(a.pos_ii, b.pos_ii, ctx.planet)
end

_terminal_rate_cap(terminal) = Inf
_terminal_rate_cap(terminal::RadioTerminal) = terminal.max_rate_bps

function link_data_rate_bps(link::DataLink, transmitter, receiver, a, b, ctx)
    rate = Float64(budget_rate_bps(link.budget, link, transmitter, receiver, norm(a.pos_ii - b.pos_ii)))
    return min(rate, _terminal_rate_cap(transmitter), _terminal_rate_cap(receiver))
end

"""
    has_data_links(model) -> Bool

Whether the interlink model `model` (or `nothing`) has a registered `DataLink`.
A run with one adds the `data_sent_bits` and `data_received_bits` states and the
data columns.
"""
has_data_links(model) = model !== nothing && any(link -> link isa DataLink, model.graph.edge_properties)

# ---- Contact windows ------------------------------------------------------------

"""
    contact_windows(history; active=false) -> Vector{NamedTuple}

The contact windows of a run's links, from its interlink history
(`model.history`): one entry `(key, start_s, end_s, duration_s)` for each span of
consecutive accepted steps over which a link was available (could close), or,
with `active=true`, selected. The history records the available links in a run
with a `DataLink`, or with `InterlinkModel(; record_available=true)`; without
them, only `active=true` finds windows. A window starts at the first step that has the
link and ends at the next step that does not, or at the last step of the
history. Windows are sorted by start time, then key.

The resolution is the scheduler's: the run's accepted steps. A resumed run takes
its history back from beside its checkpoint, retaining the earlier recorded
samples and windows. With the same configuration, checkpoint schedule and
deterministic models whose complete run state is restored, it matches the
uninterrupted checkpointed history. A run without checkpoints can take different
adaptive steps and record different
scheduling samples. Without usable
history the run warns and records only the resumed interval, where a window open
there starts; a resume that solves nothing returns an empty schedule instead.
"""
function contact_windows(history; active::Bool=false)
    return NamedTuple{(:key, :start_s, :end_s, :duration_s), Tuple{InterLinkKey, Float64, Float64, Float64}}[
        (key=window.key, start_s=window.start_s, end_s=window.end_s, duration_s=window.end_s - window.start_s)
        for window in _link_windows(history; active)]
end

# The spans of consecutive history samples over which each link was selected
# (`active=true`) or available: `(key, first, last, start_s, end_s)`, with `first`
# and `last` the first and the last sample that hold the link, `start_s` the time
# of `first`, and `end_s` the time of the next sample, or of the history's last
# sample when the link holds to the end. Sorted by start time, then key.
# `contact_windows` and the link schedule read them.
function _link_windows(history; active::Bool)
    windows = NamedTuple{(:key, :first, :last, :start_s, :end_s), Tuple{InterLinkKey, Int, Int, Float64, Float64}}[]
    isempty(history) && return windows
    open = Dict{InterLinkKey, Int}()
    for (index, sample) in enumerate(history)
        keys_now = Set{InterLinkKey}(active ? sample.active : sample.available)
        for (key, first) in collect(open)
            key in keys_now && continue
            push!(windows, (key=key, first=first, last=index - 1, start_s=history[first].time, end_s=sample.time))
            delete!(open, key)
        end
        for key in keys_now
            haskey(open, key) || (open[key] = index)
        end
    end
    finish = history[end].time
    for (key, first) in open
        push!(windows, (key=key, first=first, last=length(history), start_s=history[first].time, end_s=finish))
    end
    return sort!(windows; by=window -> (window.start_s, window.key))
end

"""
    contact_time_s(history, key; active=false) -> Float64

The total time (s) that the link `key` was available, or with `active=true`
selected, over a run's interlink history: the sum of its `contact_windows`.
"""
contact_time_s(history, key; active::Bool=false) =
    sum((window.duration_s for window in contact_windows(history; active) if window.key == key); init=0.0)
