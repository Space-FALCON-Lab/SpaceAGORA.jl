const IMPACT_ALTITUDE_M = 50_000.0

function get_impact_callback(num_sats::Int)
    function condition!(out, u, t, integrator)
        p = integrator.p
        Rp_e = p.args.environment_model.planet.Rp_e
        @inbounds for i in 1:num_sats
            out[i] = norm(_simulation_engine_module()._state_position_ii(u, i)) - Rp_e - IMPACT_ALTITUDE_M
        end
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        if p.is_active[idx]
            if callback_verbose(integrator)
                println("Impact detected for satellite $idx at time $(integrator.t) seconds at altitude <= $(IMPACT_ALTITUDE_M * 1e-3) km!")
            end
            p.is_active[idx] = false
            if _simulation_engine_module()._is_gravity_backbone_state(integrator.u)
                integrator.u.x[1].sc[idx].vel .= 0.0
            end
            if all(p.is_active .== false)
                if callback_verbose(integrator)
                    println("All satellites have impacted. Stopping simulation.")
                end
                # Terminated retcode is shared with the orbit-count stop; the cause
                # matters when reading long-mission studies, so always say which.
                println("termination_cause=impact sat=$idx t_s=$(integrator.t)")
                terminate!(integrator)
            end
        end
    end

    return VectorContinuousCallback(condition!, nothing, affect_downcrossing!, num_sats)
end

function get_orbit_end_callback(num_sats::Int)
    function condition!(out, u, t, integrator)
        # Use radial-velocity root events for orbit bookkeeping.
        # At apoapsis, dot(r,v) crosses + -> -, so -dot(r,v) crosses - -> + (upcrossing),
        # which matches the single affect! handler below.
        @inbounds for i in 1:num_sats
            pos = _simulation_engine_module()._state_position_ii(u, i)
            vel = _simulation_engine_module()._state_velocity_ii(u, i)
            out[i] = -dot(pos, vel)
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p

        p.orbit_counter[idx] += 1 # Increment the orbit counter in the shared buffers
        completed_orbits = p.orbit_counter[idx] - 1
        target_orbits = p.args.mission_configuration.number_of_orbits
        if callback_verbose(integrator)
            println("Orbit $(completed_orbits) completed by Satellite $idx at time $(integrator.t) seconds!")
        end

        # MissionOrbits should end when every active satellite reaches the requested
        # number of completed orbits, analogous to drag-passage event termination.
        if p.args.mission_configuration.mission_type == MissionOrbits && completed_orbits >= target_orbits
            all_active_reached_target = true
            @inbounds for sat_idx in eachindex(p.orbit_counter)
                if p.is_active[sat_idx] && (p.orbit_counter[sat_idx] - 1) < target_orbits
                    all_active_reached_target = false
                    break
                end
            end
            if all_active_reached_target
                if callback_verbose(integrator)
                    println("Target orbit count reached for all active satellites. Stopping simulation.")
                end
                println("termination_cause=orbit_count sat=$idx orbits=$completed_orbits t_s=$(integrator.t)")
                if applicable(terminate!, integrator)
                    terminate!(integrator)
                end
            end
        end
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end

function get_entry_end_callback(num_sats::Int, args::SimulationConfiguration)
    target_entries = _entry_target_count()
    target_entries > 0 || throw(ArgumentError("Entry target callback requires SPACEAGORA_ENTRY_TARGET_COUNT > 0"))
    entry_interface_m = args.environment_model.EI * 1e3
    entry_counter = zeros(Int64, num_sats)

    function condition!(out, u, t, integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
        @inbounds for i in 1:num_sats
            if !p.is_active[i]
                out[i] = 1.0
                continue
            end
            alt = norm(_simulation_engine_module()._state_position_ii(u, i)) - planet.Rp_e
            out[i] = alt - entry_interface_m
        end
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        if !p.is_active[idx]
            return nothing
        end

        entry_counter[idx] += 1
        completed_entries = entry_counter[idx]
        if callback_verbose(integrator)
            println("Entry $(completed_entries) detected for Satellite $idx at time $(integrator.t) seconds!")
        end

        if completed_entries >= target_entries
            all_active_reached_target = true
            @inbounds for sat_idx in eachindex(entry_counter)
                if p.is_active[sat_idx] && entry_counter[sat_idx] < target_entries
                    all_active_reached_target = false
                    break
                end
            end
            if all_active_reached_target
                if callback_verbose(integrator)
                    println("Target entry count reached for all active satellites. Stopping simulation.")
                end
                if applicable(terminate!, integrator)
                    terminate!(integrator)
                end
            end
        end
        return nothing
    end

    return VectorContinuousCallback(condition!, nothing, affect_downcrossing!, num_sats)
end

function get_drag_state_callback(num_sats::Int)
    condition!(out, u, t, integrator) = begin
        @inbounds for i in 1:num_sats
            alt = norm(_simulation_engine_module()._state_position_ii(u, i)) - integrator.p.args.environment_model.planet.Rp_e
            out[i] = alt - integrator.p.args.environment_model.EI*1e3 # Positive when above the atmosphere, negative when in the atmosphere
        end
    end
    function affect_upcrossing!(integrator, idx::Int64)
        p = integrator.p
        if callback_verbose(integrator)
            println("Switching to space integration at time $(integrator.t) seconds!")
        end
        p.shared_buffers.in_atmosphere[idx] = false
        p.shared_buffers.in_atmosphere_sample_t[idx] = Float64(integrator.t)
        # Invalidate the vacuum-predicted GRAM cache so the next atmospheric entry
        # rebuilds it from the correct state rather than interpolating stale data.
        if idx <= length(p.shared_buffers.vacuum_gram_caches)
            cache = p.shared_buffers.vacuum_gram_caches[idx]
            if cache !== nothing
                cache.valid = false
            end
        end
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_orbit # Increase the maximum timestep when exiting the atmosphere
        reltol_new, abstol_new = _callback_tolerances_for_phase(
            integrator.opts.reltol,
            integrator.opts.abstol,
            p.args,
            false
        )
        integrator.opts.reltol = reltol_new # Adjust tolerances when exiting the atmosphere
        integrator.opts.abstol = abstol_new
        schedule_event_driven_thruster_controls!(integrator, idx)
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        if callback_verbose(integrator)
            println("Switching to atmosphere integration at time $(integrator.t) seconds!")
        end
        p.shared_buffers.in_atmosphere[idx] = true
        p.shared_buffers.in_atmosphere_sample_t[idx] = Float64(integrator.t)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_atmosphere # Decrease the maximum timestep when entering the atmosphere
        reltol_new, abstol_new = _callback_tolerances_for_phase(
            integrator.opts.reltol,
            integrator.opts.abstol,
            p.args,
            true
        )
        integrator.opts.reltol = reltol_new # Adjust tolerances when entering the atmosphere
        integrator.opts.abstol = abstol_new
    end

    return VectorContinuousCallback(condition!, affect_upcrossing!, affect_downcrossing!, num_sats)
end

function get_quaternion_projection_callback(num_sats::Int, args::SimulationConfiguration)
    correction_tol = max(32 * eps(Float64), args.integration_tolerances.abstol_quaternion)
    condition(u, t, integrator) = begin
        p = integrator.p
        @inbounds for i in 1:num_sats
            if p.is_active[i]
                return true
            end
        end
        return false
    end

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        corrected = false
        @inbounds for i in 1:num_sats
            if !p.is_active[i]
                continue
            end
            q = _simulation_engine_module()._state_quaternion(u, i)
            q === nothing && continue
            qnorm2 = dot(q, q)
            if !(isfinite(qnorm2) && qnorm2 > eps(Float64)) || abs(qnorm2 - 1.0) > correction_tol
                corrected = true
            end
            # Keep accepted-step attitude states on the unit-quaternion manifold.
            u.sc[i].q .= _simulation_model_module.project_unit_quaternion(q)
        end
        if corrected && callback_verbose(integrator)
            println("Quaternion projection applied at time $(integrator.t) seconds.")
        end
    end

    # No before/after saves: the projection runs on every accepted step, and
    # the step is saved once after the discrete callbacks have run, so the
    # stored state is the projected one the next step starts from. With the
    # DiscreteCallback default the solution held the pre- and post-projection
    # states of every step as well.
    return DiscreteCallback(
        condition,
        affect!;
        initialize=(cb, u, t, integrator) -> affect!(integrator),
        save_positions=(false, false)
    )
end

function get_data_saving_callback(
    num_sats::Int,
    args::SimulationConfiguration,
    save_fields,
    saved_values=nothing
)
    saved_values = isnothing(saved_values) ? SavedValues(Float64, SaveData) : saved_values
    function save_func(u, t, integrator)
        return _save_snapshot(save_fields, u, t, integrator)
    end
    data_rate = args.mission_configuration.data_rate
    data_rate > 0.0 || throw(ArgumentError("mission_configuration.data_rate must be > 0.0, got $data_rate."))
    return SavingCallback(save_func, saved_values; saveat=data_rate, save_everystep=false)
end

"""
    TrajectoryRecorder(args; data_rate=args.mission_configuration.data_rate, save_fields=nothing)

Preallocated fixed-cadence recorder for large constellation output.

The recorder preserves the same `SaveField` payload as the existing
fixed-cadence `SavingCallback`, but stores built-in fields in typed arrays
during the solve instead of retaining one `SaveData` dictionary per sample.
Attach it with [`get_trajectory_recorder_callback`](@ref) and run with
`return_solution=false` when the caller needs fixed-cadence output but not a
full `ODESolution` history.
"""
mutable struct TrajectoryRecorder
    t::Vector{Float64}
    save_fields::Vector{SaveField}
    storage::Dict{Symbol, Any}
    fallback::Dict{Symbol, Vector{Any}}
    direct_fields::Set{Symbol}
    fused_default_fields::Bool
    count::Int
    num_sats::Int
    data_rate::Float64
end

@inline function _trajectory_recorder_capacity(args::SimulationConfiguration, data_rate::Float64)::Int
    data_rate > 0.0 || throw(ArgumentError("TrajectoryRecorder data_rate must be > 0.0, got $data_rate."))
    return max(2, ceil(Int, args.mission_configuration.mission_time / data_rate) + 2)
end

function TrajectoryRecorder(
    args::SimulationConfiguration;
    data_rate::Real=args.mission_configuration.data_rate,
    capacity::Union{Nothing, Integer}=nothing,
    save_fields=nothing,
)
    n_sats = length(args.dynamics_model.spacecraft)
    n_sats > 0 || throw(ArgumentError("TrajectoryRecorder requires at least one spacecraft."))
    rate = Float64(data_rate)
    cap = capacity === nothing ? _trajectory_recorder_capacity(args, rate) : Int(capacity)
    cap > 0 || throw(ArgumentError("TrajectoryRecorder capacity must be > 0, got $cap."))
    fields = _resolve_save_fields(save_fields, args)
    use_builtin_fillers = isnothing(save_fields)
    storage = Dict{Symbol, Any}()
    fallback = Dict{Symbol, Vector{Any}}()
    direct_fields = Set{Symbol}()
    for field in fields
        n_comp = use_builtin_fillers ? _trajectory_builtin_components(field.name) : 0
        if field.per_satellite && n_comp == 1
            storage[field.name] = Matrix{Float64}(undef, n_sats, cap)
            push!(direct_fields, field.name)
        elseif field.per_satellite && n_comp > 1
            storage[field.name] = Array{Float64, 3}(undef, n_comp, n_sats, cap)
            push!(direct_fields, field.name)
        else
            fallback[field.name] = Vector{Any}(undef, cap)
        end
    end
    return TrajectoryRecorder(
        Vector{Float64}(undef, cap),
        fields,
        storage,
        fallback,
        direct_fields,
        use_builtin_fillers,
        0,
        n_sats,
        rate,
    )
end

@inline function _trajectory_builtin_components(name::Symbol)::Int
    name === :position && return 3
    name === :velocity && return 3
    name === :wind && return 3
    name === :drag && return 3
    name === :lift && return 3
    name === :cross && return 3
    name === :quaternion && return 4
    name === :altitude && return 1
    name === :latitude_deg && return 1
    name === :longitude_deg && return 1
    name === :mass && return 1
    name === :periapsis_altitude && return 1
    name === :heat_rate && return 1
    name === :heat_load && return 1
    return 0
end

function _grow_trajectory_recorder!(rec::TrajectoryRecorder)::Nothing
    old_cap = length(rec.t)
    new_cap = max(old_cap + 1, 2 * old_cap)
    new_t = Vector{Float64}(undef, new_cap)
    if rec.count > 0
        idx = 1:rec.count
        new_t[idx] .= @view rec.t[idx]
    end
    rec.t = new_t
    for (name, data) in rec.storage
        if data isa Matrix{Float64}
            grown = Matrix{Float64}(undef, size(data, 1), new_cap)
            rec.count > 0 && (grown[:, 1:rec.count] .= @view data[:, 1:rec.count])
            rec.storage[name] = grown
        elseif data isa Array{Float64, 3}
            grown = Array{Float64, 3}(undef, size(data, 1), size(data, 2), new_cap)
            rec.count > 0 && (grown[:, :, 1:rec.count] .= @view data[:, :, 1:rec.count])
            rec.storage[name] = grown
        end
    end
    for (name, data) in rec.fallback
        grown = Vector{Any}(undef, new_cap)
        rec.count > 0 && (grown[1:rec.count] .= @view data[1:rec.count])
        rec.fallback[name] = grown
    end
    return nothing
end

"""
    reset_trajectory_recorder!(rec)

Reset `rec` so it can be reused for another solve without reallocating its
storage.
"""
function reset_trajectory_recorder!(rec::TrajectoryRecorder)::TrajectoryRecorder
    rec.count = 0
    return rec
end

@inline function _trajectory_scalar_storage(rec::TrajectoryRecorder, name::Symbol)::Matrix{Float64}
    return rec.storage[name]::Matrix{Float64}
end

@inline function _trajectory_vector_storage(rec::TrajectoryRecorder, name::Symbol)::Array{Float64, 3}
    return rec.storage[name]::Array{Float64, 3}
end

@inline function _record_vector3!(storage::Array{Float64, 3}, sat_idx::Int, sample_idx::Int, value)
    storage[1, sat_idx, sample_idx] = value[1]
    storage[2, sat_idx, sample_idx] = value[2]
    storage[3, sat_idx, sample_idx] = value[3]
    return nothing
end

function _record_default_save_fields_fused!(
    rec::TrajectoryRecorder,
    u,
    t,
    integrator,
    sample_idx::Int,
)::Nothing
    engine = _simulation_engine_module()
    p = integrator.p
    args = p.args
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    t_float = Float64(t)
    et = p.shared_buffers.et_start[] + t_float

    position_storage = _trajectory_vector_storage(rec, :position)
    velocity_storage = _trajectory_vector_storage(rec, :velocity)
    wind_storage = _trajectory_vector_storage(rec, :wind)
    drag_storage = _trajectory_vector_storage(rec, :drag)
    lift_storage = _trajectory_vector_storage(rec, :lift)
    cross_storage = _trajectory_vector_storage(rec, :cross)

    altitude_storage = _trajectory_scalar_storage(rec, :altitude)
    latitude_storage = _trajectory_scalar_storage(rec, :latitude_deg)
    longitude_storage = _trajectory_scalar_storage(rec, :longitude_deg)
    mass_storage = _trajectory_scalar_storage(rec, :mass)
    periapsis_storage = _trajectory_scalar_storage(rec, :periapsis_altitude)
    heat_rate_storage = _trajectory_scalar_storage(rec, :heat_rate)
    heat_load_storage = _trajectory_scalar_storage(rec, :heat_load)

    quaternion_storage = :quaternion in rec.direct_fields ?
        _trajectory_vector_storage(rec, :quaternion) : nothing

    winds = p.shared_buffers.winds
    drag_cache = p.save_cache.drag_cache
    lift_cache = p.save_cache.lift_cache
    cross_cache = p.save_cache.cross_cache
    shared_heat_rates = p.shared_buffers.heat_rates
    zero_vec = SVector{3, Float64}(0.0, 0.0, 0.0)

    @inbounds for sat_idx in 1:rec.num_sats
        pos = engine._state_position_ii(u, sat_idx)
        vel = engine._state_velocity_ii(u, sat_idx)
        _record_vector3!(position_storage, sat_idx, sample_idx, pos)
        _record_vector3!(velocity_storage, sat_idx, sample_idx, vel)

        rp, _ = r_intor_p!(pos, vel, planet, et, ephemerides_model)
        latlong = rtolatlong(rp, planet, ephemerides_model)
        altitude_storage[sat_idx, sample_idx] = latlong[1]
        latitude_storage[sat_idx, sample_idx] = rad2deg(latlong[2])
        longitude_storage[sat_idx, sample_idx] = rad2deg(latlong[3])

        mass_storage[sat_idx, sample_idx] = engine._state_mass_kg(u, args, sat_idx)
        _record_vector3!(wind_storage, sat_idx, sample_idx, sat_idx <= length(winds) ? winds[sat_idx] : zero_vec)
        _record_vector3!(drag_storage, sat_idx, sample_idx, sat_idx <= length(drag_cache) ? drag_cache[sat_idx] : zero_vec)
        _record_vector3!(lift_storage, sat_idx, sample_idx, sat_idx <= length(lift_cache) ? lift_cache[sat_idx] : zero_vec)
        _record_vector3!(cross_storage, sat_idx, sample_idx, sat_idx <= length(cross_cache) ? cross_cache[sat_idx] : zero_vec)

        oe = rvtoorbitalelement(pos, vel, planet)
        periapsis_storage[sat_idx, sample_idx] = oe[1] * (1.0 - oe[2]) - planet.Rp_e

        rates = if hasproperty(u, :sc)
            _compute_stage_heat_rates!(p, u.sc[sat_idx], sat_idx, t_float; use_buffered_density=false)
        elseif sat_idx <= length(shared_heat_rates)
            shared_heat_rates[sat_idx]
        else
            nothing
        end
        heat_rate_storage[sat_idx, sample_idx] =
            rates !== nothing && !isempty(rates) ? maximum(rates) : 0.0

        loads = engine._state_heat_loads(u, args, sat_idx)
        heat_load_storage[sat_idx, sample_idx] = isempty(loads) ? 0.0 : maximum(loads)

        if quaternion_storage !== nothing
            q = engine._state_quaternion(u, sat_idx)
            q === nothing && throw(ArgumentError("Quaternion save field requires orientation state."))
            quaternion_storage[1, sat_idx, sample_idx] = q[1]
            quaternion_storage[2, sat_idx, sample_idx] = q[2]
            quaternion_storage[3, sat_idx, sample_idx] = q[3]
            quaternion_storage[4, sat_idx, sample_idx] = q[4]
        end
    end

    return nothing
end

"""
    record_trajectory_sample!(rec, u, t, integrator)

Append one sample from the solver state `u` into `rec`.
"""
function record_trajectory_sample!(rec::TrajectoryRecorder, u, t, integrator)::Nothing
    rec.count == length(rec.t) && _grow_trajectory_recorder!(rec)
    sample_idx = rec.count + 1
    rec.t[sample_idx] = Float64(t)
    if rec.fused_default_fields
        _record_default_save_fields_fused!(rec, u, t, integrator, sample_idx)
    else
        for field in rec.save_fields
            rec.fallback[field.name][sample_idx] = field.getter(u, t, integrator)
        end
    end
    rec.count = sample_idx
    return nothing
end

"""
    get_trajectory_recorder_callback(rec)

Return a `SavingCallback` that records into `rec` at `rec.data_rate` seconds.
Use this as an `extra_callbacks` entry in `run_simulation`.
"""
function get_trajectory_recorder_callback(rec::TrajectoryRecorder)
    saved_values = SavedValues(Float64, Nothing)
    save_func = (u, t, integrator) -> begin
        record_trajectory_sample!(rec, u, t, integrator)
        return nothing
    end
    return SavingCallback(save_func, saved_values; saveat=rec.data_rate, save_everystep=false)
end

"""Return a view of the recorded sample times."""
trajectory_times(rec::TrajectoryRecorder) = @view rec.t[1:rec.count]

"""Return a view of one preallocated field by `SaveField.name`."""
function trajectory_field(rec::TrajectoryRecorder, name::Symbol)
    if haskey(rec.storage, name)
        data = rec.storage[name]
        if data isa Matrix{Float64}
            return @view data[:, 1:rec.count]
        elseif data isa Array{Float64, 3}
            return @view data[:, :, 1:rec.count]
        end
    elseif haskey(rec.fallback, name)
        return @view rec.fallback[name][1:rec.count]
    end
    throw(KeyError(name))
end

"""Return a `(3, num_sats, samples)` view of recorded inertial positions."""
trajectory_positions(rec::TrajectoryRecorder) = trajectory_field(rec, :position)

"""Return a `(3, num_sats, samples)` view of recorded inertial velocities."""
trajectory_velocities(rec::TrajectoryRecorder) = trajectory_field(rec, :velocity)

"""Return a `(num_sats, samples)` view of recorded spacecraft masses."""
trajectory_masses(rec::TrajectoryRecorder) = trajectory_field(rec, :mass)

function _trajectory_saved_field(rec::TrajectoryRecorder, field::SaveField, sample_idx::Int)
    if haskey(rec.storage, field.name)
        data = rec.storage[field.name]
        if data isa Matrix{Float64}
            values = Vector{Float64}(undef, rec.num_sats)
            @inbounds for sat_idx in 1:rec.num_sats
                values[sat_idx] = data[sat_idx, sample_idx]
            end
            return values
        elseif data isa Array{Float64, 3}
            n_comp = size(data, 1)
            if n_comp == 3
                values = Vector{SVector{3, Float64}}(undef, rec.num_sats)
                @inbounds for sat_idx in 1:rec.num_sats
                    values[sat_idx] = SVector{3, Float64}(
                        data[1, sat_idx, sample_idx],
                        data[2, sat_idx, sample_idx],
                        data[3, sat_idx, sample_idx],
                    )
                end
                return values
            elseif n_comp == 4
                values = Vector{SVector{4, Float64}}(undef, rec.num_sats)
                @inbounds for sat_idx in 1:rec.num_sats
                    values[sat_idx] = SVector{4, Float64}(
                        data[1, sat_idx, sample_idx],
                        data[2, sat_idx, sample_idx],
                        data[3, sat_idx, sample_idx],
                        data[4, sat_idx, sample_idx],
                    )
                end
                return values
            end
        end
    end
    return rec.fallback[field.name][sample_idx]
end

"""
    trajectory_save_data(rec)

Materialize the recorded samples as the same `Vector{SaveData}` shape produced
by the existing fixed-cadence saver. This conversion allocates at the output
boundary; it is not used during the solve.
"""
function trajectory_save_data(rec::TrajectoryRecorder)::Vector{SaveData}
    snapshots = Vector{SaveData}(undef, rec.count)
    @inbounds for sample_idx in 1:rec.count
        snapshot = SaveData()
        for field in rec.save_fields
            snapshot[field.name] = _trajectory_saved_field(rec, field, sample_idx)
        end
        snapshots[sample_idx] = snapshot
    end
    return snapshots
end
