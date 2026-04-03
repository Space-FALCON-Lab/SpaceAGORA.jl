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
        positions[i] = SVector{3, Float64}(u.sc[i].pos)
    end
    return positions
end

@inline function _save_velocities(num_sats::Int, u, t, integrator)
    velocities = Vector{SVector{3, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        velocities[i] = SVector{3, Float64}(u.sc[i].vel)
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

@inline function _save_periapsis_altitude(num_sats::Int, u, t, integrator)
    planet = integrator.p.args.environment_model.planet
    periapsis_altitudes = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        # State is in J2000; convert to PCI (body-equatorial) for rvtoorbitalelement
        pos_pci = planet.J2000_to_pci * SVector{3, Float64}(u.sc[i].pos)
        vel_pci = planet.J2000_to_pci * SVector{3, Float64}(u.sc[i].vel)
        oe = rvtoorbitalelement(pos_pci, vel_pci, planet)
        periapsis_altitudes[i] = oe[1] * (1.0 - oe[2]) - planet.Rp_e
    end
    return periapsis_altitudes
end

@inline function _save_heat_rate(num_sats::Int, u, t, integrator)
    heat_rates = Vector{Float64}(undef, num_sats)
    shared_heat_rates = integrator.p.shared_buffers.heat_rates
    @inbounds for i in 1:num_sats
        heat_rates[i] = i <= length(shared_heat_rates) ? sum(shared_heat_rates[i]) : 0.0
    end
    return heat_rates
end

@inline function _save_heat_load(num_sats::Int, u, t, integrator)
    heat_loads = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        heat_loads[i] = sum(u.sc[i].heat_loads)
    end
    return heat_loads
end

@inline function _save_mass(num_sats::Int, u, t, integrator)
    masses = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        masses[i] = Float64(u.sc[i].mass)
    end
    return masses
end

@inline function _save_quaternion(num_sats::Int, u, t, integrator)
    quaternions = Vector{SVector{4, Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        quaternions[i] = SVector{4, Float64}(u.sc[i].q)
    end
    return quaternions
end

function default_save_fields(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    fields = SaveField[
        SaveField(:position, (u, t, integrator) -> _save_positions(num_sats, u, t, integrator); per_satellite=true, column_prefix="pos"),
        SaveField(:velocity, (u, t, integrator) -> _save_velocities(num_sats, u, t, integrator); per_satellite=true, column_prefix="vel"),
        SaveField(:mass, (u, t, integrator) -> _save_mass(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:drag, (u, t, integrator) -> _save_drag(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:periapsis_altitude, (u, t, integrator) -> _save_periapsis_altitude(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:heat_rate, (u, t, integrator) -> _save_heat_rate(num_sats, u, t, integrator); per_satellite=true),
        SaveField(:heat_load, (u, t, integrator) -> _save_heat_load(num_sats, u, t, integrator); per_satellite=true)
    ]
    if args.mission_configuration.orientation_sim
        push!(fields, SaveField(:quaternion, (u, t, integrator) -> _save_quaternion(num_sats, u, t, integrator); per_satellite=true, column_prefix="q"))
    end
    return fields
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
