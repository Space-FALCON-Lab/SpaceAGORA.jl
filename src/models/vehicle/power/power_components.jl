# The power-system records: solar arrays, a battery and constant loads,
# mounted on spacecraft links with `add_power!` (`PowerSystems`). Included by
# `Components` (src/models/vehicle/spacecraft/components.jl), which owns and
# exports the records, so serialized spacecraft keep their type identity.

"""
    SolarArray(; normal, area_m2, efficiency, max_power_w=Inf)

A solar array on a spacecraft link: cells of `area_m2` square metres facing
`normal` (in the link frame, so an array on a hinged panel turns with it),
converting the incident sunlight at `efficiency`. It generates
`area_m2 * efficiency * S * max(0, cos θ) * ν` watts, capped at
`max_power_w`, where `S` is the solar flux at the spacecraft's distance from
the Sun (1361 W/m² at 1 AU), `θ` the angle between the normal and the Sun and
`ν` the lit fraction of the Sun's disk (the eclipse model of the solar
radiation pressure). Only the cell side generates.
"""
struct SolarArray
    normal::SVector{3, Float64}
    area_m2::Float64
    efficiency::Float64
    max_power_w::Float64

    function SolarArray(normal, area_m2::Real, efficiency::Real, max_power_w::Real)
        n = SVector{3, Float64}(normal)
        nn = norm(n)
        (isfinite(nn) && nn > 0.0) || throw(ArgumentError("SolarArray normal must be finite and nonzero, got $(n)."))
        (isfinite(area_m2) && area_m2 > 0.0) || throw(ArgumentError("SolarArray area_m2 must be finite and positive, got $(area_m2)."))
        (0.0 < efficiency <= 1.0) || throw(ArgumentError("SolarArray efficiency must be in (0, 1], got $(efficiency)."))
        (max_power_w > 0.0 && !isnan(max_power_w)) || throw(ArgumentError("SolarArray max_power_w must be positive, got $(max_power_w)."))
        return new(n / nn, Float64(area_m2), Float64(efficiency), Float64(max_power_w))
    end
end

SolarArray(; normal, area_m2::Real, efficiency::Real, max_power_w::Real=Inf) =
    SolarArray(normal, area_m2, efficiency, max_power_w)

"""
    Battery(; capacity_wh, initial_soc=NaN, charge_efficiency=1.0, discharge_efficiency=1.0,
            max_charge_w=Inf, max_discharge_w=Inf, recovery_soc=0.05)

The spacecraft's battery, mounted on the root link. Its stored energy, in
watt-hours, is a state the engine integrates: it charges with the surplus of
the generated power over the loads (times `charge_efficiency`, at most
`max_charge_w`) and discharges to cover a shortfall (divided by
`discharge_efficiency`, at most `max_discharge_w`). A full battery curtails
the surplus; an empty one sheds actuator load until it has recharged to
`recovery_soc` (see `PowerSystems`). `initial_soc`, a fraction of
`capacity_wh`, defaults to the spacecraft's `battery_energy_index` (a
percentage).
"""
mutable struct Battery
    capacity_wh::Float64
    initial_soc::Float64
    charge_efficiency::Float64
    discharge_efficiency::Float64
    max_charge_w::Float64
    max_discharge_w::Float64
    recovery_soc::Float64
    # Run state, switched only between integration steps (see PowerSystems):
    # full (surplus curtailed), empty (shortfall shed until recovery_soc), the
    # shedding tier, and a scratch copy of the actuator levels with the shed
    # ones zeroed.
    full::Bool
    empty::Bool
    shed_tier::Int
    shed_levels::Vector{Float64}

    function Battery(capacity_wh::Real, initial_soc::Real, charge_efficiency::Real, discharge_efficiency::Real,
                     max_charge_w::Real, max_discharge_w::Real, recovery_soc::Real)
        (isfinite(capacity_wh) && capacity_wh > 0.0) || throw(ArgumentError("Battery capacity_wh must be finite and positive, got $(capacity_wh)."))
        (isnan(initial_soc) || 0.0 <= initial_soc <= 1.0) || throw(ArgumentError("Battery initial_soc must be in [0, 1] (or NaN for the spacecraft's battery_energy_index), got $(initial_soc)."))
        (0.0 < charge_efficiency <= 1.0) || throw(ArgumentError("Battery charge_efficiency must be in (0, 1], got $(charge_efficiency)."))
        (0.0 < discharge_efficiency <= 1.0) || throw(ArgumentError("Battery discharge_efficiency must be in (0, 1], got $(discharge_efficiency)."))
        (max_charge_w > 0.0 && !isnan(max_charge_w)) || throw(ArgumentError("Battery max_charge_w must be positive, got $(max_charge_w)."))
        (max_discharge_w > 0.0 && !isnan(max_discharge_w)) || throw(ArgumentError("Battery max_discharge_w must be positive, got $(max_discharge_w)."))
        (0.0 < recovery_soc <= 1.0) || throw(ArgumentError("Battery recovery_soc must be in (0, 1], got $(recovery_soc)."))
        return new(Float64(capacity_wh), Float64(initial_soc), Float64(charge_efficiency), Float64(discharge_efficiency),
                   Float64(max_charge_w), Float64(max_discharge_w), Float64(recovery_soc), false, false, 0, Float64[])
    end
end

Battery(; capacity_wh::Real, initial_soc::Real=NaN, charge_efficiency::Real=1.0, discharge_efficiency::Real=1.0,
        max_charge_w::Real=Inf, max_discharge_w::Real=Inf, recovery_soc::Real=0.05) =
    Battery(capacity_wh, initial_soc, charge_efficiency, discharge_efficiency, max_charge_w, max_discharge_w, recovery_soc)

"""
    PowerLoad(; power_w, name="")

A constant electrical load on a spacecraft link (avionics, heaters, a
payload), in watts. Loads are never shed.
"""
struct PowerLoad
    power_w::Float64
    name::String

    function PowerLoad(power_w::Real, name::AbstractString)
        (isfinite(power_w) && power_w >= 0.0) || throw(ArgumentError("PowerLoad power_w must be finite and nonnegative, got $(power_w)."))
        return new(Float64(power_w), String(name))
    end
end

PowerLoad(; power_w::Real, name::AbstractString="") = PowerLoad(power_w, name)
