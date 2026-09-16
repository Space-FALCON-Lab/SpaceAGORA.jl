"""
    SaveField(name, getter; per_satellite=false, column_prefix=String(name))

One saved output column group. `getter(u, t, integrator)` is called at every
saved sample and returns the value: one entry per spacecraft when
`per_satellite` is true (written as `sc{i}_{column_prefix}`, with `_1.._n`
appended when the entry is a vector), or a single value for the whole run
otherwise. `name` identifies the field and must be unique within one run's save
set.

Build the built-in ones with [`save_field`](@ref) or
`default_save_fields(args; extra=...)`; construct a `SaveField` directly for a
quantity SpaceAGORA does not already know how to write.
"""
struct SaveField{F}
    name::Symbol
    getter::F
    per_satellite::Bool
    column_prefix::String
end

function SaveField(
    name::Symbol,
    getter::F;
    per_satellite::Bool=false,
    column_prefix::AbstractString=String(name)
) where {F}
    return SaveField{F}(name, getter, per_satellite, String(column_prefix))
end

@inline function _save_positions(num_sats::Int, u, t, integrator)
    positions = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        positions[i] = _simulation_engine_module()._state_position_ii(u, i)
    end
    return positions
end

@inline function _save_velocities(num_sats::Int, u, t, integrator)
    velocities = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        velocities[i] = _simulation_engine_module()._state_velocity_ii(u, i)
    end
    return velocities
end

@inline function _save_drag(num_sats::Int, u, t, integrator)
    p = integrator.p
    drag_cache = p.save_cache.drag_cache
    drags = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        drags[i] = i <= length(drag_cache) ? drag_cache[i] : SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return drags
end

@inline function _save_lift(num_sats::Int, u, t, integrator)
    p = integrator.p
    lift_cache = p.save_cache.lift_cache
    lifts = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        lifts[i] = i <= length(lift_cache) ? lift_cache[i] : SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return lifts
end

@inline function _save_cross(num_sats::Int, u, t, integrator)
    p = integrator.p
    cross_cache = p.save_cache.cross_cache
    crosses = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        crosses[i] = i <= length(cross_cache) ? cross_cache[i] : SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return crosses
end

@inline function _save_wind(num_sats::Int, u, t, integrator)
    p = integrator.p
    shared_winds = p.shared_buffers.winds
    winds = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        winds[i] = i <= length(shared_winds) ? shared_winds[i] : SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return winds
end

# Osculating classical orbital elements per satellite, recomputed from the same
# inertial state the row's position and velocity columns carry: semimajor axis
# (m), eccentricity, and the four angles in degrees. Opt-in through
# `orbital_elements_save_field`.
@inline function _save_orbital_elements(num_sats::Int, u, t, integrator)
    planet = integrator.p.args.environment_model.planet
    engine = _simulation_engine_module()
    elements = Vector{SVector{6, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        oe = rvtoorbitalelement(engine._state_position_ii(u, i), engine._state_velocity_ii(u, i), planet)
        elements[i] = SVector{6, Float64}(
            oe[1],
            oe[2],
            rad2deg(oe[3]),
            rad2deg(oe[4]),
            rad2deg(oe[5]),
            rad2deg(oe[6])
        )
    end
    return elements
end

# Inertial gravitational acceleration per satellite (m/s^2), summed over the
# run's gravity effectors only. Opt-in through `gravity_accel_save_field`.
#
# Re-evaluated here rather than cached from the right-hand side. The right-hand
# side runs every solver stage of every step, accepted or not, so caching would
# pay on the order of 10^5 stores to serve 10^3 saved samples -- and the cached
# value would be the one from the last stage the solver happened to take, at a
# different time and state than the row it landed in. One gravity evaluation per
# sample is both cheaper overall and the value that belongs in the row.
@inline function _save_gravity_accel(num_sats::Int, u, t, integrator)
    engine = _simulation_engine_module()
    p = integrator.p
    accelerations = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        accelerations[i] = engine.gravity_acceleration_ii(u, p, i, Float64(t))
    end
    return accelerations
end

@inline function _save_periapsis_altitude(num_sats::Int, u, t, integrator)
    planet = integrator.p.args.environment_model.planet
    periapsis_altitudes = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        pos = _simulation_engine_module()._state_position_ii(u, i)
        vel = _simulation_engine_module()._state_velocity_ii(u, i)
        oe = rvtoorbitalelement(pos, vel, planet)
        periapsis_altitudes[i] = oe[1] * (1.0 - oe[2]) - planet.Rp_e
    end
    return periapsis_altitudes
end

@inline function _save_altitude(num_sats::Int, u, t, integrator)
    planet = integrator.p.args.environment_model.planet
    et = integrator.p.shared_buffers.et_start[] + Float64(t)
    ephemerides_model = integrator.p.args.environment_model.ephemerides_model
    altitudes = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        pos = _simulation_engine_module()._state_position_ii(u, i)
        vel = _simulation_engine_module()._state_velocity_ii(u, i)
        rp, _ = r_intor_p!(pos, vel, planet, et, ephemerides_model)
        altitudes[i] = rtolatlong(rp, planet, ephemerides_model)[1]
    end
    return altitudes
end

@inline function _save_latitude_deg(num_sats::Int, u, t, integrator)
    planet = integrator.p.args.environment_model.planet
    et = integrator.p.shared_buffers.et_start[] + Float64(t)
    ephemerides_model = integrator.p.args.environment_model.ephemerides_model
    latitudes_deg = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        pos = _simulation_engine_module()._state_position_ii(u, i)
        vel = _simulation_engine_module()._state_velocity_ii(u, i)
        rp, _ = r_intor_p!(pos, vel, planet, et, ephemerides_model)
        latitudes_deg[i] = rad2deg(rtolatlong(rp, planet, ephemerides_model)[2])
    end
    return latitudes_deg
end

@inline function _save_longitude_deg(num_sats::Int, u, t, integrator)
    planet = integrator.p.args.environment_model.planet
    et = integrator.p.shared_buffers.et_start[] + Float64(t)
    ephemerides_model = integrator.p.args.environment_model.ephemerides_model
    longitudes_deg = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        pos = _simulation_engine_module()._state_position_ii(u, i)
        vel = _simulation_engine_module()._state_velocity_ii(u, i)
        rp, _ = r_intor_p!(pos, vel, planet, et, ephemerides_model)
        longitudes_deg[i] = rad2deg(rtolatlong(rp, planet, ephemerides_model)[3])
    end
    return longitudes_deg
end

@inline function _save_heat_rate(num_sats::Int, u, t, integrator)
    heat_rates = Vector{Float64}(undef, num_sats)
    shared_heat_rates = integrator.p.shared_buffers.heat_rates
    @inbounds for i in 1:num_sats
        rates = if hasproperty(u, :sc)
            _compute_stage_heat_rates!(integrator.p, u.sc[i], i, Float64(t); use_buffered_density=false)
        elseif i <= length(shared_heat_rates)
            shared_heat_rates[i]
        else
            Float64[]
        end
        heat_rates[i] = !isempty(rates) ? maximum(rates) : 0.0
    end
    return heat_rates
end

@inline function _save_heat_load(num_sats::Int, u, t, integrator)
    heat_loads = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        loads = _simulation_engine_module()._state_heat_loads(u, integrator.p.args, i)
        heat_loads[i] = isempty(loads) ? 0.0 : maximum(loads)
    end
    return heat_loads
end

@inline function _save_mass(num_sats::Int, u, t, integrator)
    masses = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        masses[i] = _simulation_engine_module()._state_mass_kg(u, integrator.p.args, i)
    end
    return masses
end

@inline function _save_quaternion(num_sats::Int, u, t, integrator)
    quaternions = Vector{SVector{4, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        q = _simulation_engine_module()._state_quaternion(u, i)
        q === nothing && throw(ArgumentError("Quaternion save field requires orientation state."))
        quaternions[i] = q
    end
    return quaternions
end

"""
    default_save_fields(args; extra=()) -> Vector{SaveField}

The columns every run writes, plus whatever `extra` asks for. `extra` takes
built-in field names (`Symbol`s -- see [`available_save_fields`](@ref)) and
`SaveField`s of your own, in any mix, and a name the defaults already cover is
skipped rather than duplicated:

```julia
run_simulation(args; save_fields=default_save_fields(args; extra=(:orbital_elements, :gravity_accel)))
```
"""
function default_save_fields(args::SimulationConfiguration; extra=())
    num_sats = length(args.dynamics_model.spacecraft)
    fields = SaveField[
        SaveField(:position, (u, t, integrator) -> _save_positions(num_sats, u, t, integrator); per_satellite=true, column_prefix="pos"),
        SaveField(:velocity, (u, t, integrator) -> _save_velocities(num_sats, u, t, integrator); per_satellite=true, column_prefix="vel"),
        SaveField(:altitude, (u, t, integrator) -> _save_altitude(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:latitude_deg, (u, t, integrator) -> _save_latitude_deg(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:longitude_deg, (u, t, integrator) -> _save_longitude_deg(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:mass, (u, t, integrator) -> _save_mass(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:wind, (u, t, integrator) -> _save_wind(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:drag, (u, t, integrator) -> _save_drag(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:lift, (u, t, integrator) -> _save_lift(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:cross, (u, t, integrator) -> _save_cross(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:periapsis_altitude, (u, t, integrator) -> _save_periapsis_altitude(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:heat_rate, (u, t, integrator) -> _save_heat_rate(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:heat_load, (u, t, integrator) -> _save_heat_load(num_sats, u, t, integrator); per_satellite=true)
    ]
    if args.mission_configuration.orientation_sim
        push!(fields, quaternion_save_field(args))
    end
    for requested in extra
        field = requested isa SaveField ? requested : save_field(Symbol(requested), args)
        any(existing -> existing.name === field.name, fields) && continue
        push!(fields, field)
    end
    return fields
end

"""
    orbital_elements_save_field(args) -> SaveField

The `orbital_elements` field: each spacecraft's osculating classical elements,
written as `sc{i}_orbital_elements_1..6` in the order semimajor axis (m),
eccentricity, inclination (deg), right ascension of the ascending node (deg),
argument of periapsis (deg), true anomaly (deg).
"""
function orbital_elements_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField(:orbital_elements, (u, t, integrator) -> _save_orbital_elements(num_sats, u, t, integrator); per_satellite=true, column_prefix="orbital_elements")
end

"""
    gravity_accel_save_field(args) -> SaveField

The `gravity_accel` field: the inertial gravitational acceleration on each
spacecraft (m/s^2), written as `sc{i}_gravity_accel_1..3`. Summed over the
run's gravity effectors only -- constant, inverse-square, J2, spherical
harmonics and third-body -- so drag, SRP and thrust are excluded.
"""
function gravity_accel_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField(:gravity_accel, (u, t, integrator) -> _save_gravity_accel(num_sats, u, t, integrator); per_satellite=true, column_prefix="gravity_accel")
end

"""
    quaternion_save_field(args) -> SaveField

The `quaternion` field: each spacecraft's inertial-to-body attitude quaternion,
written as `sc{i}_q_1..4` (scalar-last). Part of `default_save_fields` whenever
the run integrates attitude; requesting it from a run with
`orientation_sim = false` throws at the first saved sample, because there is no
attitude state to read.
"""
function quaternion_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField(:quaternion, (u, t, integrator) -> _save_quaternion(num_sats, u, t, integrator); per_satellite=true, column_prefix="q")
end

# Built-in fields addressable by name. Each entry builds its field from the
# run's configuration, so a name is resolved in exactly one place no matter
# whether it arrives through `default_save_fields(; extra=...)` or
# `save_field`.
const _SAVE_FIELD_BUILDERS = Dict{Symbol, Function}(
    :orbital_elements => orbital_elements_save_field,
    :gravity_accel => gravity_accel_save_field,
    :quaternion => quaternion_save_field,
)

"""
    available_save_fields() -> Vector{Symbol}

Every built-in save field that can be requested by name, sorted. Pass any of
them to `default_save_fields(args; extra=...)` or [`save_field`](@ref).
"""
available_save_fields()::Vector{Symbol} = sort!(collect(keys(_SAVE_FIELD_BUILDERS)))

"""
    save_field(name::Symbol, args) -> SaveField

Build one built-in save field by name, for the run `args` describes. Throws on
an unknown name, listing what is available.
"""
function save_field(name::Symbol, args::SimulationConfiguration)::SaveField
    builder = get(_SAVE_FIELD_BUILDERS, name, nothing)
    builder === nothing && throw(ArgumentError(
        "Unknown save field $(repr(name)). Built-in names: $(join(available_save_fields(), ", ")). " *
        "Pass a SaveField of your own for anything else."
    ))
    return builder(args)
end

@inline function _resolve_save_fields(save_fields, args::SimulationConfiguration)
    resolved = isnothing(save_fields) ? default_save_fields(args) : collect(save_fields)
    names = Symbol[field.name for field in resolved]
    length(unique(names)) == length(names) || throw(ArgumentError("save_fields names must be unique. Got $(names)."))
    return resolved
end

function _save_snapshot(save_fields, u, t, integrator)::SaveData
    # SaveData is a persistence/output boundary; keep runtime logic on typed state and buffers.
    snapshot = SaveData()
    for field in save_fields
        snapshot[field.name] = field.getter(u, t, integrator)
    end
    return snapshot
end
