@inline function _evaluate_dynamic_effector(
    effector,
    sc_view,
    state_sample,
    p,
    sat_idx::Int,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if _wrench_method_available(effector)
        state_sample === nothing && error("wrench-based effector evaluation requires a typed StateSample.")
        req = SimulationModel.environment_requirements(effector)
        env = if req.atmosphere && p.shared_buffers.rhs_atmosphere_prefilled[]
            sample_environment_with_buffered_atm(req, effector, sc_view, p, sat_idx, t)
        else
            sample_environment(req, effector, sc_view, p, sat_idx, t; write_buffers=false)
        end
        return SimulationModel.wrench(effector, state_sample, env, t)
    end
    return SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)
end

@inline function _solver_partition_validated(effector)::Symbol
    partition = SimulationModel.solver_partition(effector)
    if partition === :implicit || partition === :explicit
        return partition
    end
    throw(ArgumentError(
        "solver_partition($(nameof(typeof(effector)))) must return :implicit or :explicit, got $(repr(partition))."
    ))
end

@inline _effector_in_partition(effector, partition::Symbol)::Bool =
    _solver_partition_validated(effector) === partition

@inline function _partition_needs_state_sample(
    dynamic_effectors::Tuple,
    partition::Symbol,
)::Bool
    return any(_effector_in_partition(effector, partition) && _wrench_method_available(effector) for effector in dynamic_effectors)
end

@inline function _partition_selected_count(
    dynamic_effectors::Tuple,
    partition::Symbol,
)::Int
    count = 0
    @inbounds for effector in dynamic_effectors
        count += _effector_in_partition(effector, partition) ? 1 : 0
    end
    return count
end

@inline function _accumulate_dynamic_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    dynamic_effectors::Tuple,
    effector_decision
)
    # Only pay the time_ns() syscall overhead when telemetry is actually needed.
    needs_timing = effector_decision.policy_applied
    effector_started_ns = needs_timing ? time_ns() : UInt64(0)
    n_effectors = length(dynamic_effectors)
    needs_state_sample = any(_wrench_method_available(effector) for effector in dynamic_effectors)
    state_sample = if needs_state_sample
        spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
        build_state_sample(sc_view, spacecraft, p.args.mission_configuration.orientation_sim)
    else
        nothing
    end
    if effector_decision.use_threads
        reduced = SimulationModel.ParallelPolicy.threaded_reduce(
            n_effectors,
            effector_decision.allotment,
            () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (local_sum, eff_idx) -> begin
                effector = dynamic_effectors[eff_idx]
                force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
                local_sum[1] += force[1]
                local_sum[2] += force[2]
                local_sum[3] += force[3]
                local_sum[4] += torque[1]
                local_sum[5] += torque[2]
                local_sum[6] += torque[3]
                return nothing
            end,
            (dest, src) -> begin
                @inbounds for i in 1:6
                    dest[i] += src[i]
                end
                return nothing
            end
        )
        # Unpack reduced MVector directly into forces/torques without an intermediate SVector.
        @inbounds begin
            forces[1] = reduced[1]
            forces[2] = reduced[2]
            forces[3] = reduced[3]
            torques[1] = reduced[4]
            torques[2] = reduced[5]
            torques[3] = reduced[6]
        end
    else
        @inbounds for effector in dynamic_effectors
            force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
            forces .+= force
            torques .+= torque
        end
    end
    if needs_timing
        elapsed_ns = Int64(time_ns() - effector_started_ns)
        if sat_idx == 1
            _update_effector_cost_model!(p.shared_buffers, n_effectors, elapsed_ns, effector_decision.allotment)
        end
        SimulationModel.ParallelPolicy.record_policy_observation!(
            :dynamic_effectors;
            mode=effector_decision.mode,
            num_items=n_effectors,
            use_threads=effector_decision.use_threads,
            elapsed_ns=elapsed_ns
        )
    end
    return nothing
end

@inline function _accumulate_dynamic_effectors_partitioned!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    dynamic_effectors::Tuple,
    effector_decision,
    partition::Symbol,
)
    selected_count = _partition_selected_count(dynamic_effectors, partition)
    selected_count == 0 && return nothing

    state_sample = if _partition_needs_state_sample(dynamic_effectors, partition)
        spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
        build_state_sample(sc_view, spacecraft, p.args.mission_configuration.orientation_sim)
    else
        nothing
    end

    if effector_decision.use_threads && selected_count > 1
        reduced = SimulationModel.ParallelPolicy.threaded_reduce(
            length(dynamic_effectors),
            effector_decision.allotment,
            () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (local_sum, eff_idx) -> begin
                effector = dynamic_effectors[eff_idx]
                _effector_in_partition(effector, partition) || return nothing
                force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
                local_sum[1] += force[1]
                local_sum[2] += force[2]
                local_sum[3] += force[3]
                local_sum[4] += torque[1]
                local_sum[5] += torque[2]
                local_sum[6] += torque[3]
                return nothing
            end,
            (dest, src) -> begin
                @inbounds for i in 1:6
                    dest[i] += src[i]
                end
                return nothing
            end
        )
        forces .= SVector{3, Float64}(reduced[1], reduced[2], reduced[3])
        torques .= SVector{3, Float64}(reduced[4], reduced[5], reduced[6])
        return nothing
    end

    @inbounds for effector in dynamic_effectors
        _effector_in_partition(effector, partition) || continue
        force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
        forces .+= force
        torques .+= torque
    end
    return nothing
end

@inline function _accumulate_dynamic_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    dynamic_effectors::Tuple,
    effector_decision
)
    t = if p === nothing
        0.0
    else
        Float64(p.shared_buffers.current_time[])
    end
    return _accumulate_dynamic_effectors!(forces, torques, sc_view, p, sat_idx, t, dynamic_effectors, effector_decision)
end

@inline function _accumulate_control_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    debug_control::Bool
)::Float64
    mass_rate = 0.0
    @inbounds for control_effector in p.args.control_model.control_effectors
        control_force, control_torque = SimulationModel.calcControlForceTorque(control_effector, sc_view, p, sat_idx, t)
        control_mass_rate = SimulationModel.calcControlMassFlowRate(control_effector, sc_view, p, sat_idx, t)
        if debug_control && (norm(control_force) > 0.0 || norm(control_torque) > 0.0)
            println("Applying control effect for spacecraft $sat_idx at time $t seconds:")
            println("  Control force: $control_force")
        end
        forces .+= control_force
        torques .+= control_torque
        mass_rate += isfinite(control_mass_rate) ? control_mass_rate : 0.0
    end
    return mass_rate
end

@inline function _ensure_rhs_flat_effector_scratch!(
    shared_buffers,
    num_sats::Int,
    workers::Int,
)
    partials_ref = shared_buffers.rhs_flat_effector_partials
    partials = partials_ref[]
    if size(partials, 1) != 6 || size(partials, 2) < num_sats || size(partials, 3) < workers
        partials_ref[] = zeros(Float64, 6, num_sats, workers)
    else
        partials[:, 1:num_sats, 1:workers] .= 0.0
    end

    totals_ref = shared_buffers.rhs_flat_effector_totals
    totals = totals_ref[]
    if size(totals, 1) != 8 || size(totals, 2) < num_sats
        totals_ref[] = zeros(Float64, 8, num_sats)
    else
        totals[:, 1:num_sats] .= 0.0
    end

    samples_ref = shared_buffers.rhs_flat_state_samples
    samples = samples_ref[]
    if length(samples) < num_sats
        resize!(samples, num_sats)
    end
    @inbounds for i in 1:num_sats
        samples[i] = nothing
    end

    items_ref = shared_buffers.rhs_flat_work_items
    items = items_ref[]
    if length(items) < num_sats
        resize!(items, num_sats)
    end
    return nothing
end

@inline function _flat_partition_selected(
    effector,
    partition::Union{Nothing, Symbol},
)::Bool
    partition === nothing && return true
    return _effector_in_partition(effector, partition)
end

@inline function _rhs_effector_cost_rank(effector)::Int
    if effector isa SimulationModel.GravitationalHarmonicsModel
        return 4
    elseif effector isa SimulationModel.NBodyGravityModel
        return 4
    elseif effector isa SimulationModel.SolarRadiationPressureModel
        return 3
    elseif effector isa SimulationModel.AerodynamicCoefficientfM
        return 2
    elseif effector isa SimulationModel.InverseSquaredJ2GravityModel
        return 2
    elseif effector isa SimulationModel.InverseSquaredGravityModel
        return 1
    end
    return 1
end

@inline function _rhs_flat_item_cost_key(dynamic_effectors::Tuple, item::Int)::Int
    n_effectors = length(dynamic_effectors)
    eff_idx = mod(item - 1, n_effectors) + 1
    return -_rhs_effector_cost_rank(dynamic_effectors[eff_idx])
end

function _prepare_rhs_flat_work_items!(
    work_items::Vector{Int},
    p,
    dynamic_effectors::Tuple,
    num_sats::Int,
    partition::Union{Nothing, Symbol},
)::Int
    n_effectors = length(dynamic_effectors)
    required = num_sats * n_effectors
    if length(work_items) < required
        resize!(work_items, required)
    end

    count_items = 0
    @inbounds for sat_idx in 1:num_sats
        p.is_active[sat_idx] || continue
        base = (sat_idx - 1) * n_effectors
        for eff_idx in 1:n_effectors
            effector = dynamic_effectors[eff_idx]
            _flat_partition_selected(effector, partition) || continue
            count_items += 1
            work_items[count_items] = base + eff_idx
        end
    end
    if count_items > 1
        sort!(@view(work_items[1:count_items]); by=item -> _rhs_flat_item_cost_key(dynamic_effectors, item), alg=Base.Sort.DEFAULT_STABLE)
    end
    return count_items
end

function _accumulate_harmonics_flat_batch!(
    sc_state,
    p,
    t::Float64,
    model::SimulationModel.GravitationalHarmonicsModel,
    plan,
)::Nothing
    num_sats = length(sc_state)
    active_sats = count(identity, p.is_active)
    active_sats <= 0 && return nothing
    workers = SimulationModel.ParallelPolicy.thread_worker_count(active_sats, plan.allotment)
    _ensure_rhs_flat_effector_scratch!(p.shared_buffers, num_sats, workers)
    totals = p.shared_buffers.rhs_flat_effector_totals[]
    work_items = p.shared_buffers.rhs_flat_work_items[]
    if length(work_items) < active_sats
        resize!(work_items, active_sats)
    end

    count_items = 0
    @inbounds for sat_idx in 1:num_sats
        p.is_active[sat_idx] || continue
        count_items += 1
        work_items[count_items] = sat_idx
    end

    et = p.shared_buffers.et_start[] + t
    lpi = SimulationModel.DynamicEffectors.PerturbationEffectors._harmonics_lpi_at!(model, p, et)

    needs_timing = plan.policy_applied
    started_ns = needs_timing ? time_ns() : UInt64(0)

    # Batched SIMD kernel: partition work_items into n_workers contiguous ranges.
    # Each worker runs _harmonics_flat_batch_kernel! over its slice, loading
    # coefficients once per (degree,order) pair and iterating over the batch
    # with @turbo SIMD over the satellite dimension.
    n_workers = min(workers, count_items)
    batch_size = cld(count_items, max(1, n_workers))
    pool = SimulationModel.DynamicEffectors.PerturbationEffectors._get_harmonics_batch_pool(model, n_workers, batch_size)
    if n_workers <= 1
        SimulationModel.DynamicEffectors.PerturbationEffectors._harmonics_flat_batch_kernel!(
            totals, model, sc_state, work_items, 1, count_items, lpi, pool[1]
        )
    else
        Threads.@sync for w in 1:n_workers
            item_start = (w - 1) * batch_size + 1
            item_end   = min(w * batch_size, count_items)
            item_start > count_items && continue
            ws = pool[w]
            Threads.@spawn SimulationModel.DynamicEffectors.PerturbationEffectors._harmonics_flat_batch_kernel!(
                totals, model, sc_state, work_items, item_start, item_end, lpi, ws
            )
        end
    end

    if needs_timing
        elapsed_ns = Int64(time_ns() - started_ns)
        _update_effector_cost_model!(p.shared_buffers, max(1, active_sats), elapsed_ns, plan.allotment)
        SimulationModel.ParallelPolicy.record_policy_observation!(
            :dynamic_effectors;
            mode=:flat_constellation_effector_queue,
            num_items=max(1, active_sats),
            use_threads=true,
            elapsed_ns=elapsed_ns,
        )
    end
    return nothing
end

function _accumulate_dynamic_effectors_flat_batch!(
    sc_state,
    p,
    t::Float64,
    dynamic_effectors::Tuple,
    plan;
    partition::Union{Nothing, Symbol}=nothing,
)
    num_sats = length(sc_state)
    n_effectors = length(dynamic_effectors)
    num_items = num_sats * n_effectors
    num_items <= 0 && return nothing
    workers = SimulationModel.ParallelPolicy.thread_worker_count(num_items, plan.allotment)
    _ensure_rhs_flat_effector_scratch!(p.shared_buffers, num_sats, workers)
    selected_count = partition === nothing ? n_effectors : _partition_selected_count(dynamic_effectors, partition)
    selected_count <= 0 && return nothing

    if partition === nothing &&
       n_effectors == 1 &&
       dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel
        return _accumulate_harmonics_flat_batch!(sc_state, p, t, dynamic_effectors[1], plan)
    end

    partials = p.shared_buffers.rhs_flat_effector_partials[]
    totals = p.shared_buffers.rhs_flat_effector_totals[]
    state_samples = p.shared_buffers.rhs_flat_state_samples[]
    work_items = p.shared_buffers.rhs_flat_work_items[]

    needs_state_sample = if partition === nothing
        any(_wrench_method_available(effector) for effector in dynamic_effectors)
    else
        _partition_needs_state_sample(dynamic_effectors, partition)
    end
    if needs_state_sample
        spacecraft = p.args.dynamics_model.spacecraft
        orientation_sim = p.args.mission_configuration.orientation_sim
        @inbounds for sat_idx in 1:num_sats
            if p.is_active[sat_idx]
                @views state_samples[sat_idx] = build_state_sample(sc_state[sat_idx], spacecraft[sat_idx], orientation_sim)
            end
        end
    end

    needs_timing = plan.policy_applied
    started_ns = needs_timing ? time_ns() : UInt64(0)
    count_items = _prepare_rhs_flat_work_items!(
        work_items,
        p,
        dynamic_effectors,
        num_sats,
        partition,
    )
    count_items <= 0 && return nothing

    SimulationModel.ParallelPolicy.threaded_foreach_worker(count_items, plan.allotment) do worker_id, item_idx
        item = @inbounds work_items[item_idx]
        sat_idx = fld(item - 1, n_effectors) + 1
        eff_idx = mod(item - 1, n_effectors) + 1
        effector = dynamic_effectors[eff_idx]
        @views sc_view = sc_state[sat_idx]
        state_sample = needs_state_sample ? state_samples[sat_idx] : nothing
        force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
        @inbounds begin
            partials[1, sat_idx, worker_id] += force[1]
            partials[2, sat_idx, worker_id] += force[2]
            partials[3, sat_idx, worker_id] += force[3]
            partials[4, sat_idx, worker_id] += torque[1]
            partials[5, sat_idx, worker_id] += torque[2]
            partials[6, sat_idx, worker_id] += torque[3]
        end
        return nothing
    end

    @inbounds for sat_idx in 1:num_sats
        for worker_id in 1:workers
            totals[1, sat_idx] += partials[1, sat_idx, worker_id]
            totals[2, sat_idx] += partials[2, sat_idx, worker_id]
            totals[3, sat_idx] += partials[3, sat_idx, worker_id]
            totals[4, sat_idx] += partials[4, sat_idx, worker_id]
            totals[5, sat_idx] += partials[5, sat_idx, worker_id]
            totals[6, sat_idx] += partials[6, sat_idx, worker_id]
        end
    end

    if needs_timing
        elapsed_ns = Int64(time_ns() - started_ns)
        active_sats = count(identity, p.is_active)
        _update_effector_cost_model!(
            p.shared_buffers,
            max(1, count_items),
            elapsed_ns,
            plan.allotment,
        )
        SimulationModel.ParallelPolicy.record_policy_observation!(
            :dynamic_effectors;
            mode=:flat_constellation_effector_queue,
            num_items=max(1, count_items),
            use_threads=true,
            elapsed_ns=elapsed_ns,
        )
    end
    return nothing
end

@inline function _flat_totals_force_torque(totals, sat_idx::Int)
    forces = MVector{3, Float64}(
        totals[1, sat_idx],
        totals[2, sat_idx],
        totals[3, sat_idx],
    )
    torques = MVector{3, Float64}(
        totals[4, sat_idx],
        totals[5, sat_idx],
        totals[6, sat_idx],
    )
    return forces, torques
end

# Warm SpiceRhsMemo with solar and N-body positions for time t before the parallel region
# starts. Workers then find immediate cache hits and only hold the memo lock for a fast
# dict lookup rather than an expensive SPICE kernel call.
function _prefill_shared_body_samples!(p, t::Float64, sc_state, dynamic_effectors::Tuple)::Nothing
    first_active = findfirst(p.is_active)
    first_active === nothing && return nothing
    @views sc_view = sc_state[first_active]

    for effector in dynamic_effectors
        if effector isa SimulationModel.SolarRadiationPressureModel
            sample_solar_ephemeris(sc_view, p, first_active, t)
            break
        end
    end

    for effector in dynamic_effectors
        if effector isa SimulationModel.NBodyGravityModel
            sample_third_body_ephemerides(effector, sc_view, p, first_active, t)
            break
        end
    end
    return nothing
end

# Parallel atmosphere pre-sample: fills shared_buffers.densities/temperatures/winds for
# all active satellites before the flat effector queue. l_pi is pre-computed once so the
# batch avoids acquiring harmonics_lpi_lock per satellite.
function _prefill_atmosphere_samples!(p, t::Float64, sc_state)::Nothing
    num_sats = length(sc_state)
    num_sats == 0 && return nothing
    planet = p.args.environment_model.planet
    l_pi = _planet_lpi_at_engine(p, t)
    @batch for sat_idx in 1:num_sats
        p.is_active[sat_idx] || continue
        @views sc_view = sc_state[sat_idx]
        planet_frame = sample_planet_frame_with_lpi(sc_view, planet, l_pi)
        _sample_atmosphere_from_planet_frame(sc_view, planet_frame, p, sat_idx, t; write_buffers=true)
    end
    return nothing
end

function _spacecraft_dynamics_flat_constellation_effector_queue!(
    du::ComponentVector,
    u::ComponentVector,
    p,
    t::Float64,
    plan;
    rhs_kind::Symbol,
    partition::Union{Nothing, Symbol}=nothing,
)
    sc_state = u.sc
    sc_du = du.sc
    spacecraft = p.args.dynamics_model.spacecraft
    dynamic_effectors = p.args.dynamics_model.dynamic_effectors
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t

    # Phase 1 (serial): warm SpiceRhsMemo for solar/N-body bodies before workers start.
    _prefill_shared_body_samples!(p, t, sc_state, dynamic_effectors)

    # Phase 2 (parallel): pre-sample atmosphere per satellite so the hot loop reads buffers.
    needs_atm_prefill = partition === nothing && any(dynamic_effectors) do eff
        _wrench_method_available(eff) && SimulationModel.environment_requirements(eff).atmosphere
    end
    if needs_atm_prefill
        _prefill_atmosphere_samples!(p, t, sc_state)
        p.shared_buffers.rhs_atmosphere_prefilled[] = true
    end

    _accumulate_dynamic_effectors_flat_batch!(sc_state, p, t, dynamic_effectors, plan; partition=partition)

    if needs_atm_prefill
        p.shared_buffers.rhs_atmosphere_prefilled[] = false
    end
    totals = p.shared_buffers.rhs_flat_effector_totals[]

    @inbounds for i in eachindex(sc_state)
        if !p.is_active[i]
            sc_du[i] .= 0.0
            continue
        end
        @views begin
            sc_view = sc_state[i]
            du_view = sc_du[i]
            forces, torques = _flat_totals_force_torque(totals, i)
            if rhs_kind == :implicit
                SimulationModel.DynamicsTranslational.assign_force_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )
                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end
                du_view.heat_loads .= 0.0
            elseif rhs_kind == :slow
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=true,
                )
                SimulationModel.DynamicsTranslational.assign_slow_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )
                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end
                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            else
                mass_rate = rhs_kind == :explicit || rhs_kind == :full ?
                    _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control) :
                    0.0
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=true,
                )
                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )
                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end
                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
    return nothing
end

@inline function _assign_heat_rate_derivative!(du_heat::AbstractVector, heat_rates::AbstractVector)
    if length(heat_rates) == length(du_heat)
        du_heat .= heat_rates
        return nothing
    end
    du_heat .= 0.0
    n_copy = min(length(heat_rates), length(du_heat))
    @inbounds for j in 1:n_copy
        du_heat[j] = heat_rates[j]
    end
    return nothing
end

@inline function _gravity_backbone_xyz_chunk(state, sat_idx::Int)::SVector{3, Float64}
    if hasproperty(state, :sc)
        sc_state = state.sc[sat_idx]
        if hasproperty(sc_state, :pos)
            return SVector{3, Float64}(sc_state.pos)
        elseif hasproperty(sc_state, :vel)
            return SVector{3, Float64}(sc_state.vel)
        end
        return SVector{3, Float64}(sc_state)
    end

    if sat_idx <= length(state)
        sat_state = state[sat_idx]
        if hasproperty(sat_state, :pos)
            return SVector{3, Float64}(sat_state.pos)
        elseif hasproperty(sat_state, :vel)
            return SVector{3, Float64}(sat_state.vel)
        end
    end

    base = 3 * (sat_idx - 1)
    return SVector{3, Float64}(state[base + 1], state[base + 2], state[base + 3])
end

@inline function _gravity_backbone_state_sample(q_state, dq_state, p, sat_idx::Int)::StateSample
    spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
    pos_ii = _gravity_backbone_xyz_chunk(q_state, sat_idx)
    vel_ii = _gravity_backbone_xyz_chunk(dq_state, sat_idx)
    mass_kg = spacecraft.dry_mass + spacecraft.prop_mass
    return StateSample(pos_ii, vel_ii, mass_kg; spacecraft=spacecraft)
end

@inline function _gravity_backbone_state_sample(u_state, p, sat_idx::Int)::StateSample
    return _gravity_backbone_state_sample(
        _gravity_backbone_position_state(u_state),
        _gravity_backbone_velocity_state(u_state),
        p,
        sat_idx,
    )
end

@inline function _gravity_backbone_core_acceleration(
    dynamic_effectors::Tuple,
    state_sample::StateSample,
    p,
    sat_idx::Int,
    t::Float64,
)::SVector{3, Float64}
    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for effector in dynamic_effectors
        _gravity_backbone_structure_validated(effector) === :position_only_static_gravity || continue
        req = SimulationModel.environment_requirements(effector)
        env = sample_environment(req, effector, state_sample, p, sat_idx, t; write_buffers=false)
        accel_ii .+= SimulationModel.gravity_backbone_acceleration_ii(effector, state_sample, env, t)
    end
    return SVector{3, Float64}(accel_ii)
end

@inline function _gravity_backbone_kick_acceleration(
    dynamic_effectors::Tuple,
    state_sample::StateSample,
    p,
    sat_idx::Int,
    t::Float64,
)::SVector{3, Float64}
    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for effector in dynamic_effectors
        _gravity_backbone_kick_structure_validated(effector) === :velocity_kick_explicit || continue
        req = SimulationModel.environment_requirements(effector)
        env = sample_environment(req, effector, state_sample, p, sat_idx, t; write_buffers=false)
        accel_ii .+= SimulationModel.gravity_backbone_kick_acceleration_ii(effector, state_sample, env, t)
    end
    return SVector{3, Float64}(accel_ii)
end

function _gravity_backbone_half_kick!(u_state, p, t::Float64, half_dt::Float64)
    half_dt > 0.0 || return nothing
    dynamic_effectors = p.args.dynamics_model.dynamic_effectors
    vel_state = _gravity_backbone_velocity_state(u_state)
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(vel_state.sc))

    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(vel_state.sc) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(vel_state.sc)
            if !p.is_active[i]
                continue
            end
            state_sample = _gravity_backbone_state_sample(u_state, p, i)
            accel_ii = _gravity_backbone_kick_acceleration(dynamic_effectors, state_sample, p, i, t)
            vel_state.sc[i].vel .+= half_dt .* accel_ii
        end
    else
        @inbounds for i in eachindex(vel_state.sc)
            if !p.is_active[i]
                continue
            end
            state_sample = _gravity_backbone_state_sample(u_state, p, i)
            accel_ii = _gravity_backbone_kick_acceleration(dynamic_effectors, state_sample, p, i, t)
            vel_state.sc[i].vel .+= half_dt .* accel_ii
        end
    end
    return nothing
end

function spacecraft_dynamics_gravity_backbone!(ddu, dq, q, p, t::Float64)
    q_state = q.sc
    dq_state = dq.sc
    ddu_state = ddu.sc
    dynamic_effectors = p.args.dynamics_model.dynamic_effectors
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(q_state))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(q_state) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(q_state)
            if !p.is_active[i]
                ddu_state[i].vel .= 0.0
                continue
            end
            state_sample = _gravity_backbone_state_sample(q, dq, p, i)
            ddu_state[i].vel .= _gravity_backbone_core_acceleration(dynamic_effectors, state_sample, p, i, t)
        end
    else
        @inbounds for i in eachindex(q_state)
            if !p.is_active[i]
                ddu_state[i].vel .= 0.0
                continue
            end
            state_sample = _gravity_backbone_state_sample(q, dq, p, i)
            ddu_state[i].vel .= _gravity_backbone_core_acceleration(dynamic_effectors, state_sample, p, i, t)
        end
    end
    return nothing
end

@inline function _assign_orientation_rhs!(
    du_view,
    sc_view,
    inertia_tensor::AbstractMatrix{<:Real},
    torques::AbstractVector{<:Real};
    propagate_quaternion::Bool,
    include_gyroscopic::Bool,
)
    omega_body = SimulationModel.DynamicsRotational.body_angular_velocity(sc_view.ω)
    tau_body = SimulationModel.DynamicsRotational.body_torque(torques)

    if propagate_quaternion
        du_view.q .= SimulationModel.DynamicsRotational.quaternion_derivative(omega_body, sc_view.q)
    else
        du_view.q .= 0.0
    end

    du_view.ω .= SimulationModel.DynamicsRotational.angular_acceleration(
        omega_body,
        inertia_tensor,
        tau_body;
        include_gyroscopic=include_gyroscopic,
    )
    return nothing
end

function spacecraft_dynamics!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    plan = _rhs_execution_plan(p.args, p, dynamic_effectors, length(spacecraft))
    if plan.mode == :flat_constellation_effector_queue
        return _spacecraft_dynamics_flat_constellation_effector_queue!(du, u, p, t, plan; rhs_kind=:full)
    end
    effector_decision = plan.effector_decision
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics!

function spacecraft_dynamics_slow!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    p.shared_buffers.current_time[] = t
    plan = _rhs_execution_plan(p.args, p, dynamic_effectors, length(spacecraft))
    if plan.mode == :flat_constellation_effector_queue
        return _spacecraft_dynamics_flat_constellation_effector_queue!(du, u, p, t, plan; rhs_kind=:slow)
    end
    effector_decision = plan.effector_decision
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_slow_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_slow_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics_slow!

function spacecraft_dynamics_implicit_atmosphere!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    p.shared_buffers.current_time[] = t
    plan = _rhs_execution_plan(p.args, p, dynamic_effectors, length(spacecraft))
    if plan.mode == :flat_constellation_effector_queue
        return _spacecraft_dynamics_flat_constellation_effector_queue!(
            du,
            u,
            p,
            t,
            plan;
            rhs_kind=:implicit,
            partition=:implicit,
        )
    end
    effector_decision = plan.effector_decision
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :implicit)

                SimulationModel.DynamicsTranslational.assign_force_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :implicit)

                SimulationModel.DynamicsTranslational.assign_force_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    end
end # function spacecraft_dynamics_implicit_atmosphere!

function spacecraft_dynamics_explicit_remainder!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    plan = _rhs_execution_plan(p.args, p, dynamic_effectors, length(spacecraft))
    if plan.mode == :flat_constellation_effector_queue
        return _spacecraft_dynamics_flat_constellation_effector_queue!(
            du,
            u,
            p,
            t,
            plan;
            rhs_kind=:explicit,
            partition=:explicit,
        )
    end
    effector_decision = plan.effector_decision
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :explicit)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :explicit)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics_explicit_remainder!

function spacecraft_dynamics_fast_control!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    spacecraft = p.args.dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                SimulationModel.DynamicsTranslational.assign_control_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                SimulationModel.DynamicsTranslational.assign_control_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    end
end # function spacecraft_dynamics_fast_control!

function build_initial_conditions(args)::ComponentVector
    # Each spacecraft can have a different body count, so the state layout has
    # to be sized per spacecraft before we pack the global ComponentVector.
    sc_shapes = map(args.dynamics_model.spacecraft) do sc
        n_bodies = length(sc.links)
        mass = sc.dry_mass + sc.prop_mass
        if args.mission_configuration.orientation_sim
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies),
                q = Float64[0.0, 0.0, 0.0, 1.0], 
                ω = zeros(3)
            )
        else
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies)
            )
        end
    end

    state = ComponentVector(sc = sc_shapes)

    for i in eachindex(args.dynamics_model.spacecraft)
        spacecraft = args.dynamics_model.spacecraft[i]
        sc_view = state.sc[i]
        ic = spacecraft.initial_condition
        r0, v0 = if ic isa CartesianInitialCondition
            collect(ic.pos), collect(ic.vel)   # assumed J2000 frame
        else
            # Orbit-element initial conditions are propagated directly in J2000.
            r_j2000, v_j2000 = orbitalelemtorv(ic, args.environment_model.planet)
            collect(r_j2000), collect(v_j2000)
        end
        sc_view.pos .= r0
        sc_view.vel .= v0
        sc_view.heat_loads .= 0.0  
        
        if args.mission_configuration.orientation_sim
            sc_view.q .= SimulationModel.project_unit_quaternion(spacecraft.initial_condition.q)
            sc_view.ω .= spacecraft.initial_condition.ang_vel
        end
    end

    return state
end
