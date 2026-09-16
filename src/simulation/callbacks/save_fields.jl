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

# Non-root link poses for the visualization sidecar. These are not integrated
# state on the standard path; control code rotates the `Link` objects in place,
# so the snapshot has to read them at save time. Ragged across spacecraft
# (7 floats per non-root link), so the output layer falls back to its generic
# per-satellite column expansion: `sc{i}_link_pose_{1..7n}`.
@inline function _save_link_poses(num_sats::Int, u, t, integrator)
    spacecraft = integrator.p.args.dynamics_model.spacecraft
    poses = Vector{Vector{Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        poses[i] = link_pose_vector(spacecraft[i])
    end
    return poses
end

# Atmospheric density the RHS last evaluated for each satellite (kg/m^3), for
# the viewer's pass coloring; zero outside the atmosphere or without one.
@inline function _save_density(num_sats::Int, u, t, integrator)
    densities = integrator.p.shared_buffers.densities
    out = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        out[i] = i <= length(densities) ? Float64(densities[i]) : 0.0
    end
    return out
end

# Unit vector from the planet's center to the Sun in the inertial frame of the
# saved positions, for the viewer's sun lighting. One direction per row, not
# per satellite: every spacecraft of a run shares the central body.
@inline function _save_sun_direction(u, t, integrator)
    environment = integrator.p.args.environment_model
    et = integrator.p.shared_buffers.et_start[] + Float64(t)
    direction = ephemerides_sun_direction_ii(environment.planet, et, environment.ephemerides_model)
    return direction === nothing ? SVector{3, Float64}(NaN, NaN, NaN) : direction
end

# Unit vector from the planet's center to Earth, for the viewer's earthshine.
# One direction per row like the Sun's, and only written away from Earth.
@inline function _save_earth_direction(u, t, integrator)
    environment = integrator.p.args.environment_model
    et = integrator.p.shared_buffers.et_start[] + Float64(t)
    direction = ephemerides_body_direction_ii(environment.planet, "earth", et, environment.ephemerides_model)
    return direction === nothing ? SVector{3, Float64}(NaN, NaN, NaN) : direction
end

# Cloth robot-arm chain: the integrated arm state (`arm_r`, `arm_q` per arm
# link) relative to the spacecraft position, 7 floats per link.
@inline function _save_arm_poses(num_sats::Int, u, t, integrator)
    poses = Vector{Vector{Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        sc_view = hasproperty(u, :sc) ? u.sc[i] : nothing
        r_ii = _simulation_engine_module()._state_position_ii(u, i)
        poses[i] = sc_view === nothing ? Float64[] : arm_pose_vector(sc_view, r_ii)
    end
    return poses
end

# Thruster firing levels (0 to 1) for the viewer's plumes, in the order the
# visualization scene lists the thrusters: the spacecraft's links in order,
# each link's `thrusters` in order. Ragged across spacecraft (one value per
# thruster), so the output layer expands them as `sc{i}_thruster_level_{k}`.
@inline function _thruster_count(model)::Int
    n = 0
    for link in model.links
        n += length(link.thrusters)
    end
    return n
end

@inline function _save_thruster_levels(num_sats::Int, counts::Vector{Int}, u, t, integrator)
    effectors = integrator.p.args.control_model.control_effectors
    out = Vector{Vector{Float64}}(undef, num_sats)
    @inbounds for i in 1:num_sats
        levels = zeros(Float64, counts[i])
        if counts[i] > 0
            for effector in effectors
                reported = control_thruster_levels(effector, i)
                reported === nothing && continue
                for k in 1:min(counts[i], length(reported))
                    value = Float64(reported[k])
                    levels[k] = isfinite(value) ? clamp(value, 0.0, 1.0) : 0.0
                end
                break   # the first effector that drives this spacecraft's thrusters owns them
            end
        end
        out[i] = levels
    end
    return out
end

"""
    thruster_level_counts(args) -> Vector{Int}

Number of thrusters per spacecraft, in the scene's order. Zero for a
spacecraft whose links carry none or whose control effectors report no firing
levels, so the `thruster_level` field covers only the vehicles that fly them.
"""
function thruster_level_counts(args::SimulationConfiguration)::Vector{Int}
    spacecraft = args.dynamics_model.spacecraft
    effectors = hasproperty(args, :control_model) && hasproperty(args.control_model, :control_effectors) ?
        args.control_model.control_effectors : ()
    counts = zeros(Int, length(spacecraft))
    for i in eachindex(spacecraft)
        n = _thruster_count(spacecraft[i])
        n == 0 && continue
        any(effector -> control_thruster_levels(effector, i) !== nothing, effectors) || continue
        counts[i] = n
    end
    return counts
end

"""
    thruster_level_save_field(args) -> SaveField

The `thruster_level` field: every thruster's firing level (0 to 1) per
spacecraft, written to `sc{i}_thruster_level_{k}`.
"""
function thruster_level_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    counts = thruster_level_counts(args)
    return SaveField(:thruster_level, (u, t, integrator) -> _save_thruster_levels(num_sats, counts, u, t, integrator); per_satellite=true, column_prefix="thruster_level")
end

@inline _thruster_level_field_enabled(args::SimulationConfiguration)::Bool = any(>(0), thruster_level_counts(args))

@inline function _arm_pose_field_enabled(args::SimulationConfiguration)::Bool
    args.simulation_settings.save_visualization_scene || return false
    return any(i -> robot_arm_plan_for(args, i) !== nothing, eachindex(args.dynamics_model.spacecraft))
end

# The Sun direction is available when the run's ephemerides backend can
# actually resolve the Sun for this body at the start epoch: SPICE with the
# planet's kernels furnished, or the simple model at Earth. The probe is the
# resolution itself, so a SPICE-backed configuration whose kernels were never
# loaded omits the columns instead of throwing at the first save.
function _sun_direction_field_enabled(args::SimulationConfiguration)::Bool
    return try
        environment = args.environment_model
        et = ephemerides_time_seconds(args.initial_time, environment.ephemerides_model)
        direction = ephemerides_sun_direction_ii(environment.planet, et, environment.ephemerides_model)
        direction !== nothing && all(isfinite, direction)
    catch
        false
    end
end

# Earth's direction is written only away from Earth (at Earth it is the origin
# of the frame, not a light source) and only when the ephemerides resolve it,
# probed the same way the Sun's is.
function _earth_direction_field_enabled(args::SimulationConfiguration)::Bool
    return try
        environment = args.environment_model
        lowercase(strip(String(environment.planet.name))) == "earth" && return false
        et = ephemerides_time_seconds(args.initial_time, environment.ephemerides_model)
        direction = ephemerides_body_direction_ii(environment.planet, "earth", et, environment.ephemerides_model)
        direction !== nothing && all(isfinite, direction)
    catch
        false
    end
end

@inline function _density_field_enabled(args::SimulationConfiguration)::Bool
    args.simulation_settings.save_visualization_scene || return false
    return !(args.environment_model.density_model isa NoAtmosphereModel)
end

@inline function _link_pose_field_enabled(args::SimulationConfiguration)::Bool
    args.simulation_settings.save_visualization_scene || return false
    return any(model -> !isempty(link_pose_link_indices(model)), args.dynamics_model.spacecraft)
end

# Plume-surface interaction: the descent engine's footprint on the regolith,
# read straight out of the effector's per-spacecraft state. Present only when a
# `PlumeSurfaceInteractionModel` is among the run's dynamic effectors.
@inline function _plume_effector(args::SimulationConfiguration)
    for effector in args.dynamics_model.dynamic_effectors
        effector isa PlumeSurfaceInteractionModel && return effector
    end
    return nothing
end

@inline function _save_plume(model, field::Symbol, num_sats::Int)
    values = getfield(model.state, field)
    out = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        out[i] = i <= length(values) ? Float64(values[i]) : 0.0
    end
    return out
end

"""
    plume_save_fields(args) -> Vector{SaveField}

The seven `sc{i}_plume_*` columns (height above the ground along the engine
axis, peak surface pressure and wall shear stress, mass erosion rate, eroded
mass, ejecta speed and ground-effect force) when the run carries a
`PlumeSurfaceInteractionModel`, and no columns otherwise.
"""
function plume_save_fields(args::SimulationConfiguration)
    model = _plume_effector(args)
    model === nothing && return SaveField[]
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField[
        SaveField(:plume_height_m, (u, t, integrator) -> _save_plume(model, :height_m, num_sats); per_satellite=true, column_prefix="plume_height_m"),
        SaveField(:plume_shear_pa, (u, t, integrator) -> _save_plume(model, :shear_pa, num_sats); per_satellite=true, column_prefix="plume_shear_pa"),
        SaveField(:plume_pressure_pa, (u, t, integrator) -> _save_plume(model, :pressure_pa, num_sats); per_satellite=true, column_prefix="plume_pressure_pa"),
        SaveField(:plume_erosion_kg_s, (u, t, integrator) -> _save_plume(model, :erosion_kg_s, num_sats); per_satellite=true, column_prefix="plume_erosion_kg_s"),
        SaveField(:plume_eroded_kg, (u, t, integrator) -> _save_plume(model, :eroded_kg, num_sats); per_satellite=true, column_prefix="plume_eroded_kg"),
        SaveField(:plume_ejecta_mps, (u, t, integrator) -> _save_plume(model, :ejecta_mps, num_sats); per_satellite=true, column_prefix="plume_ejecta_mps"),
        SaveField(:plume_ground_effect_n, (u, t, integrator) -> _save_plume(model, :ground_effect_n, num_sats); per_satellite=true, column_prefix="plume_ground_effect_n"),
    ]
end

function default_save_fields(args::SimulationConfiguration)
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
        push!(fields, SaveField(:quaternion, (u, t, integrator) -> _save_quaternion(num_sats, u, t, integrator); per_satellite=true, column_prefix="q"))
    end
    for field in plume_save_fields(args)
        push!(fields, field)
    end
    _thruster_level_field_enabled(args) && push!(fields, thruster_level_save_field(args))
    if _sun_direction_field_enabled(args)
        push!(fields, sun_direction_save_field(args))
    end
    if _earth_direction_field_enabled(args)
        push!(fields, earth_direction_save_field(args))
    end
    for field in visualization_save_fields(args)
        push!(fields, field)
    end
    return fields
end

"""
    density_save_field(args) -> SaveField

The `density` field (kg/m^3 per satellite) the viewer colors passes by.
"""
function density_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField(:density, (u, t, integrator) -> _save_density(num_sats, u, t, integrator); per_satellite=true, column_prefix="density")
end

"""
    sun_direction_save_field(args) -> SaveField

The `sun_dir` field: the unit vector from the planet's center to the Sun in
the inertial frame of the saved positions, written as `sun_dir_1..3` (one
direction per row, shared by every spacecraft). The viewer's sun lighting
reads it; `_sun_direction_field_enabled` decides whether the run's ephemerides
can supply it.
"""
function sun_direction_save_field(::SimulationConfiguration)
    return SaveField(:sun_dir, _save_sun_direction; per_satellite=false, column_prefix="sun_dir")
end

"""
    earth_direction_save_field(args) -> SaveField

The `earth_dir` field: the unit vector from the planet's center to Earth in the
inertial frame of the saved positions, written as `earth_dir_1..3` (one
direction per row, shared by every spacecraft). The viewer's earthshine light
reads it; `_earth_direction_field_enabled` decides whether the run's planet and
ephemerides can supply it.
"""
function earth_direction_save_field(::SimulationConfiguration)
    return SaveField(:earth_dir, _save_earth_direction; per_satellite=false, column_prefix="earth_dir")
end

"""
    visualization_save_fields(args) -> Vector{SaveField}

The extra fields the viewer wants when `save_visualization_scene` is on:
`link_pose` for articulated spacecraft and `density` when the run has an
atmosphere. Empty when the flag is off. The engine appends any of these that
an explicit `save_fields` list lacks.
"""
function visualization_save_fields(args::SimulationConfiguration)
    fields = SaveField[]
    _link_pose_field_enabled(args) && push!(fields, link_pose_save_field(args))
    _arm_pose_field_enabled(args) && push!(fields, arm_pose_save_field(args))
    _density_field_enabled(args) && push!(fields, density_save_field(args))
    return fields
end

"""
    arm_pose_save_field(args) -> SaveField

The `arm_pose` field: per spacecraft, every cloth robot-arm link's COM
position relative to the spacecraft (inertial, meters) and inertial
quaternion.
"""
function arm_pose_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField(:arm_pose, (u, t, integrator) -> _save_arm_poses(num_sats, u, t, integrator); per_satellite=true, column_prefix="arm_pose")
end

"""
    link_pose_save_field(args) -> SaveField

The `link_pose` field on its own, so a caller that passes explicit
`save_fields` to `run_simulation` still gets it when the visualization flag
is on (the engine appends it if absent).
"""
function link_pose_save_field(args::SimulationConfiguration)
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField(:link_pose, (u, t, integrator) -> _save_link_poses(num_sats, u, t, integrator); per_satellite=true, column_prefix="link_pose")
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
