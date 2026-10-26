"""
    PowerSystems

The spacecraft power budget: solar arrays generate, the battery stores, and
constant loads and the mounted actuators draw. Parts are mounted on links with
[`add_power!`](@ref). When a battery is mounted, the engine integrates its
stored energy and, when it is empty, sheds actuator load in tiers: laser
terminals first, then torque rods, then the reaction wheels' motors.
"""
module PowerSystems

using LinearAlgebra
using StaticArrays
using ..Components: SolarArray, Battery, PowerLoad, ReactionWheel, Magnetorquer, LaserTerminal, PowerReceiver
using ..SpacecraftModels: Link, SpacecraftModel
using ..ActuatorMounting: actuator_context, applied_actuator_levels, mounted_actuators, mounted_thruster_count,
    mounted_wheel_count, mounted_magnetorquer_count, mounted_laser_count
using ..Actuators: actuator_consumption
using ..Components: wheel_power

export add_power!, has_battery, has_power_parts, battery, initial_battery_wh
export solar_generation_w, bus_load_w, actuator_power_by_tier, shed_tier, battery_rate_wh_s
export shed_levels, power_balance, battery_boundary!, update_power_state!, mean_solar_distance_m
export SOLAR_FLUX_1AU_W_M2, AU_M, solar_flux_w_m2

"""Solar flux at 1 astronomical unit, W/m²."""
const SOLAR_FLUX_1AU_W_M2 = 1361.0
"""One astronomical unit, m."""
const AU_M = 1.495978707e11

"""Solar flux at `distance_m` from the Sun, W/m²."""
@inline solar_flux_w_m2(distance_m::Real)::Float64 = SOLAR_FLUX_1AU_W_M2 * (AU_M / Float64(distance_m))^2

_check_own_link(model::SpacecraftModel, link::Link) = any(l -> l === link, model.links) || throw(ArgumentError(
    "Power parts mount on one of the spacecraft's own links; this link is not in its link list."))

"""
    add_power!(model::SpacecraftModel, link::Link, part) -> part

Mount a power-system part on `link`: a `SolarArray` (any link), the `Battery`
(the root link, at most one) or a `PowerLoad` (any link). Returns the part.
"""
function add_power!(model::SpacecraftModel, link::Link, array::SolarArray)
    _check_own_link(model, link)
    push!(link.solar_arrays, array)
    return array
end

function add_power!(model::SpacecraftModel, link::Link, b::Battery)
    _check_own_link(model, link)
    link === model.root || throw(ArgumentError("The battery mounts on the spacecraft's root link."))
    isempty(link.batteries) || throw(ArgumentError("The spacecraft already has a battery; mount at most one."))
    push!(link.batteries, b)
    return b
end

function add_power!(model::SpacecraftModel, link::Link, load::PowerLoad)
    _check_own_link(model, link)
    push!(link.power_loads, load)
    return load
end

"""Whether `model` carries a battery, whose stored energy the engine integrates."""
has_battery(model::SpacecraftModel)::Bool = !isempty(model.root.batteries)

"""Whether `model` carries any power-system part: a battery, a solar array, a load, or a `PowerReceiver`."""
has_power_parts(model::SpacecraftModel)::Bool = has_battery(model) ||
    any(link -> !isempty(link.solar_arrays) || !isempty(link.power_loads) ||
        any(actuator -> actuator isa PowerReceiver, link.custom_actuators), model.links)

"""The spacecraft's battery, or `nothing`."""
battery(model::SpacecraftModel) = has_battery(model) ? model.root.batteries[1] : nothing

"""
    initial_battery_wh(model) -> Float64

The battery's stored energy at the start of a run: its capacity times its
`initial_soc`, or, when that is NaN, times the spacecraft's
`battery_energy_index` over 100.
"""
function initial_battery_wh(model::SpacecraftModel)::Float64
    b = battery(model)
    b === nothing && return 0.0
    soc = isnan(b.initial_soc) ? model.battery_energy_index / 100.0 : b.initial_soc
    return b.capacity_wh * soc
end

"""
    solar_generation_w(model, sun_dir_body, flux_w_m2, lit_fraction) -> Float64

The power the solar arrays of `model` generate with the Sun along the unit
vector `sun_dir_body` (body frame, from the spacecraft), a solar flux of
`flux_w_m2` and a lit fraction `lit_fraction` of the Sun's disk.
"""
function solar_generation_w(model::SpacecraftModel, sun_dir_body::SVector{3, Float64}, flux_w_m2::Float64, lit_fraction::Float64)::Float64
    lit_fraction > 0.0 || return 0.0
    total = 0.0
    for link in model.links
        isempty(link.solar_arrays) && continue
        to_body = actuator_context(link).link_to_body
        for array in link.solar_arrays
            c = dot(to_body * array.normal, sun_dir_body)
            c > 0.0 || continue
            total += min(array.max_power_w, array.area_m2 * array.efficiency * flux_w_m2 * c * lit_fraction)
        end
    end
    return total
end

"""The constant loads of `model`, W."""
function bus_load_w(model::SpacecraftModel)::Float64
    total = 0.0
    for link in model.links, load in link.power_loads
        total += load.power_w
    end
    return total
end

"""
    actuator_power_by_tier(model, levels; wheel_momentum=nothing) -> NTuple{4, Float64}

The mounted actuators' power at their applied levels, split by shedding tier:
(laser terminals, torque rods, reaction wheels, every other kind), W.
"""
function actuator_power_by_tier(model::SpacecraftModel, levels::AbstractVector; wheel_momentum=nothing)
    applied = applied_actuator_levels(model, levels; wheel_momentum=wheel_momentum)
    assembly = model.root.rw_assembly
    first_mounted = assembly.n_wheels - length(assembly.wheels)
    lasers = rods = wheels = other = 0.0
    j = 0
    for (k, entry) in enumerate(mounted_actuators(model))
        a = entry.actuator
        if a isa ReactionWheel
            j += 1
            h = wheel_momentum === nothing ? NaN : Float64(wheel_momentum[first_mounted + j])
            wheels += wheel_power(a, applied[k], h)
        else
            w = Float64(actuator_consumption(a, applied[k], actuator_context(entry.link)).power_w)
            if a isa LaserTerminal
                lasers += w
            elseif a isa Magnetorquer
                rods += w
            else
                other += w
            end
        end
    end
    return (lasers, rods, wheels, other)
end

"""
    shed_tier(available_w, bus_w, by_tier) -> Int

The fewest shedding tiers that bring the load within `available_w`: 0 sheds
nothing, 1 the laser terminals, 2 also the torque rods, 3 also the reaction
wheels' motors. Constant loads and other kinds are never shed; tier 3 is the
most that can be done.
"""
function shed_tier(available_w::Float64, bus_w::Float64, by_tier::NTuple{4, Float64})::Int
    lasers, rods, wheels, other = by_tier
    load = bus_w + lasers + rods + wheels + other
    load <= available_w && return 0
    load - lasers <= available_w && return 1
    load - lasers - rods <= available_w && return 2
    return 3
end

"""
    battery_rate_wh_s(b::Battery, net_w) -> Float64

The rate of change of the stored energy, Wh/s, for a net power `net_w`
(generation minus load, W): a surplus charges at `charge_efficiency` up to
`max_charge_w`, a shortfall discharges over `discharge_efficiency` up to
`max_discharge_w`; a full battery takes no surplus and an empty one gives
nothing.
"""
function battery_rate_wh_s(b::Battery, net_w::Float64)::Float64
    if net_w >= 0.0
        b.full && return 0.0
        return min(net_w, b.max_charge_w) * b.charge_efficiency / 3600.0
    end
    b.empty && return 0.0
    return -min(-net_w, b.max_discharge_w) / b.discharge_efficiency / 3600.0
end

"""
    shed_levels(model, levels) -> AbstractVector

`levels` (one per mounted actuator) with the actuators the battery's shedding
tier sheds set to zero: tier 1 the laser terminals, tier 2 also the torque
rods, tier 3 also the reaction wheels. `levels` itself when nothing is shed.
"""
function shed_levels(model::SpacecraftModel, levels::AbstractVector)
    b = battery(model)
    (b === nothing || b.shed_tier == 0) && return levels
    out = b.shed_levels
    resize!(out, length(levels))
    copyto!(out, levels)
    wheels_from = mounted_thruster_count(model) + 1
    rods_from = wheels_from + mounted_wheel_count(model)
    lasers_from = rods_from + mounted_magnetorquer_count(model)
    lasers_to = lasers_from + mounted_laser_count(model) - 1
    for k in lasers_from:lasers_to
        out[k] = 0.0
    end
    if b.shed_tier >= 2
        for k in rods_from:(lasers_from - 1)
            out[k] = 0.0
        end
    end
    if b.shed_tier >= 3
        for k in wheels_from:(rods_from - 1)
            out[k] = 0.0
        end
    end
    return out
end

"""
    power_balance(model, levels, generation_w; wheel_momentum=nothing)
        -> (load_w, unshed_load_w, by_tier, net_w)

The load on the bus with the current shedding (constant loads plus the
mounted actuators at `levels`, `nothing` for none), the load it would be
without shedding, the unshed actuator power by tier (see
[`actuator_power_by_tier`](@ref)), and the net power `generation_w - load_w`.
"""
function power_balance(model::SpacecraftModel, levels, generation_w::Float64; wheel_momentum=nothing)
    bus = bus_load_w(model)
    if levels === nothing
        return (load_w = bus, unshed_load_w = bus, by_tier = (0.0, 0.0, 0.0, 0.0), net_w = generation_w - bus)
    end
    unshed = actuator_power_by_tier(model, levels; wheel_momentum=wheel_momentum)
    applied = shed_levels(model, levels)
    shed = applied === levels ? unshed : actuator_power_by_tier(model, applied; wheel_momentum=wheel_momentum)
    load = bus + sum(shed)
    return (load_w = load, unshed_load_w = bus + sum(unshed), by_tier = unshed, net_w = generation_w - load)
end

# Whether `energy_wh` (Wh) stored in battery `b` counts as full: within `1e-9`
# of its capacity. It is the one test of fullness: the boundary event, the start
# of a run and the end of each step use it, and so does the energy-balancing
# scheduling policy.
within_full_tolerance(b, energy_wh::Real) = energy_wh >= b.capacity_wh * (1.0 - 1.0e-9)

"""
    battery_boundary!(model, energy_wh) -> Bool

At a battery boundary event: mark the battery empty when its stored energy
has fallen to zero, full when it has reached its capacity. Returns whether a
flag changed.
"""
function battery_boundary!(model::SpacecraftModel, energy_wh::Real)::Bool
    b = battery(model)
    b === nothing && return false
    tol = 1.0e-9 * b.capacity_wh
    if energy_wh <= tol && !b.empty
        b.empty = true
        return true
    elseif within_full_tolerance(b, energy_wh) && !b.full
        b.full = true
        return true
    end
    return false
end

"""
    update_power_state!(model, energy_wh, generation_w, by_tier) -> Bool

At the end of an integration step: an empty battery that has recharged to its
`recovery_soc` stops shedding, a full one that has discharged takes charge
again, and the shedding tier becomes the fewest tiers that keep the load
within what the arrays and the battery can supply (the arrays alone while the
battery is empty). Returns whether anything changed, so the caller marks the
derivative discontinuous.
"""
function update_power_state!(model::SpacecraftModel, energy_wh::Real, generation_w::Float64, by_tier::NTuple{4, Float64})::Bool
    b = battery(model)
    b === nothing && return false
    changed = false
    if b.empty && energy_wh >= b.recovery_soc * b.capacity_wh
        b.empty = false
        changed = true
    end
    if b.full && !within_full_tolerance(b, energy_wh)
        b.full = false
        changed = true
    end
    available = generation_w + (b.empty ? 0.0 : b.max_discharge_w)
    tier = shed_tier(available, bus_load_w(model), by_tier)
    if tier != b.shed_tier
        b.shed_tier = tier
        changed = true
    end
    return changed
end

# Mean distances from the Sun of the bodies with planet models, m (Titan and
# its parent Saturn share one).
const _MEAN_SOLAR_DISTANCE_AU = Dict("earth" => 1.0, "moon" => 1.0, "mars" => 1.5237, "venus" => 0.7233,
                                     "titan" => 9.5371, "saturn" => 9.5371, "jupiter" => 5.2029, "mercury" => 0.3871)

"""
    mean_solar_distance_m(planet) -> Float64

The central body's mean distance from the Sun, which sets the solar flux on
the arrays: 1 AU for Earth and the Moon, the mean orbital distance for the
other planets (Titan: Saturn's); 1 AU for a body not in the table.
"""
mean_solar_distance_m(planet)::Float64 = AU_M * get(_MEAN_SOLAR_DISTANCE_AU, lowercase(String(planet.name)), 1.0)

end # module PowerSystems
