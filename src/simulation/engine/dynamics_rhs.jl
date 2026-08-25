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
        env = sample_environment_with_reusable_buffers(req, effector, sc_view, p, sat_idx, t)
        return SimulationModel.wrench_caching!(effector, state_sample, env, t, p, sat_idx)
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

# Type-stable replacement for `for effector in dynamic_effectors`.
#
# `dynamic_effectors` is a heterogeneous tuple whose eltype is the abstract
# AbstractForceTorqueModel, so iterating it leaves the loop variable abstract:
# _evaluate_dynamic_effector becomes a dynamic dispatch and the force/torque it
# returns are inferred as Any, which boxes on every accumulate. Measured
# directly: 240 bytes per iteration for the 3-effector (gravity + SRP + aero)
# mix, versus 0 bytes for a homogeneous 1-tuple. That per-satellite,
# per-RHS-call cost is a large part of why the atmosphere constellation cases
# allocate ~30 GB per solve while the single-effector vacuum cases allocate
# ~1.6 GB, and allocation is what caps their thread scaling.
#
# Peeling one element at a time via Base.tail gives every recursion level a
# concretely typed `first(effs)`, so the chain unrolls at compile time and stays
# allocation-free. Effector tuples are small (1-5 entries), so recursion depth is
# not a concern.
@inline _accumulate_effector_chain!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    ::Tuple{},
    sc_view,
    state_sample,
    p,
    sat_idx::Int,
    t::Float64,
)::Nothing = nothing

@inline function _accumulate_effector_chain!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    effs::Tuple,
    sc_view,
    state_sample,
    p,
    sat_idx::Int,
    t::Float64,
)::Nothing
    force, torque = _evaluate_dynamic_effector(first(effs), sc_view, state_sample, p, sat_idx, t)
    # Broadcast, matching the loop this replaced: the type-stability win comes
    # from peeling the tuple, not from rewriting the accumulate, and keeping the
    # original expression avoids gratuitously perturbing floating-point rounding.
    forces .+= force
    torques .+= torque
    return _accumulate_effector_chain!(
        forces, torques, Base.tail(effs), sc_view, state_sample, p, sat_idx, t
    )
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
        _accumulate_effector_chain!(
            forces, torques, dynamic_effectors, sc_view, state_sample, p, sat_idx, t
        )
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
            elapsed_ns=elapsed_ns,
            env=_policy_env_config(p)
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
    rw_torque_body::MVector{3, Float64},
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
        rw_torque = SimulationModel.calcReactionWheelTorque(control_effector, sc_view, p, sat_idx, t)
        if debug_control && (norm(control_force) > 0.0 || norm(control_torque) > 0.0)
            println("Applying control effect for spacecraft $sat_idx at time $t seconds:")
            println("  Control force: $control_force")
        end
        forces .+= control_force
        torques .+= control_torque
        rw_torque === nothing || (rw_torque_body .+= rw_torque)
        mass_rate += isfinite(control_mass_rate) ? control_mass_rate : 0.0
    end
    return mass_rate
end

@inline function _ensure_rhs_flat_effector_scratch!(
    shared_buffers,
    num_sats::Int,
    workers::Int;
    zero_partials::Bool=true,
)
    if zero_partials
        partials_ref = shared_buffers.rhs_flat_effector_partials
        partials = partials_ref[]
        if size(partials, 1) != 6 || size(partials, 2) < num_sats || size(partials, 3) < workers
            partials_ref[] = zeros(Float64, 6, num_sats, workers)
        elseif size(partials, 2) == num_sats && size(partials, 3) == workers
            # Exact-size buffer (the steady-state case): contiguous fill! is a
            # straight memset, cheaper than the strided view broadcast below.
            fill!(partials, 0.0)
        else
            partials[:, 1:num_sats, 1:workers] .= 0.0
        end
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

    state_pos_ref = shared_buffers.rhs_flat_state_pos_ii
    state_vel_ref = shared_buffers.rhs_flat_state_vel_ii
    state_mass_ref = shared_buffers.rhs_flat_state_mass_kg
    state_q_ref = shared_buffers.rhs_flat_state_q_ib
    state_omega_ref = shared_buffers.rhs_flat_state_omega_body
    if length(state_pos_ref[]) < num_sats
        resize!(state_pos_ref[], num_sats)
        resize!(state_vel_ref[], num_sats)
        resize!(state_mass_ref[], num_sats)
        resize!(state_q_ref[], num_sats)
        resize!(state_omega_ref[], num_sats)
    end

    planet_pos_ref = shared_buffers.rhs_flat_planet_pos_pp
    planet_vel_ref = shared_buffers.rhs_flat_planet_vel_pp
    planet_alt_ref = shared_buffers.rhs_flat_planet_alt_m
    planet_lat_ref = shared_buffers.rhs_flat_planet_lat_rad
    planet_lon_ref = shared_buffers.rhs_flat_planet_lon_rad
    if length(planet_pos_ref[]) < num_sats
        resize!(planet_pos_ref[], num_sats)
        resize!(planet_vel_ref[], num_sats)
        resize!(planet_alt_ref[], num_sats)
        resize!(planet_lat_ref[], num_sats)
        resize!(planet_lon_ref[], num_sats)
    end

    items_ref = shared_buffers.rhs_flat_work_items
    items = items_ref[]
    if length(items) < num_sats
        resize!(items, num_sats)
    end

    packet_starts_ref = shared_buffers.rhs_flat_packet_starts
    packet_ends_ref = shared_buffers.rhs_flat_packet_ends
    packet_costs_ref = shared_buffers.rhs_flat_packet_costs
    packet_elapsed_ref = shared_buffers.rhs_flat_packet_elapsed_ns
    if length(packet_starts_ref[]) < num_sats
        resize!(packet_starts_ref[], num_sats)
        resize!(packet_ends_ref[], num_sats)
        resize!(packet_costs_ref[], num_sats)
        resize!(packet_elapsed_ref[], num_sats)
    end
    return nothing
end

function _prefill_rhs_flat_state_samples!(shared_buffers, sc_state, p)::Nothing
    orientation_sim = p.args.mission_configuration.orientation_sim
    pos = shared_buffers.rhs_flat_state_pos_ii[]
    vel = shared_buffers.rhs_flat_state_vel_ii[]
    mass = shared_buffers.rhs_flat_state_mass_kg[]
    q = shared_buffers.rhs_flat_state_q_ib[]
    omega = shared_buffers.rhs_flat_state_omega_body[]
    @inbounds for sat_idx in eachindex(sc_state)
        p.is_active[sat_idx] || continue
        @views sc_view = sc_state[sat_idx]
        pos_i, vel_i = _extract_sample_pos_vel(sc_view)
        pos[sat_idx] = pos_i
        vel[sat_idx] = vel_i
        mass[sat_idx] = _extract_sample_mass_kg(sc_view)
        if orientation_sim
            q[sat_idx] = SVector{4, Float64}(sc_view.q)
            omega[sat_idx] = SVector{3, Float64}(sc_view.ω)
        end
    end
    return nothing
end

@inline function _rhs_flat_state_sample_from_buffers(shared_buffers, spacecraft, sat_idx::Int, orientation_sim::Bool)
    return StateSample(
        shared_buffers.rhs_flat_state_pos_ii[][sat_idx],
        shared_buffers.rhs_flat_state_vel_ii[][sat_idx],
        shared_buffers.rhs_flat_state_mass_kg[][sat_idx];
        q_ib=orientation_sim ? shared_buffers.rhs_flat_state_q_ib[][sat_idx] : nothing,
        ω_body=orientation_sim ? shared_buffers.rhs_flat_state_omega_body[][sat_idx] : nothing,
        spacecraft=spacecraft[sat_idx],
    )
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

# Takes the per-item cost default as a value rather than reading it from the
# environment. The env reader (`_effector_cost_ns_per_item_default`) goes through
# _parse_positive_float_env, which allocates unconditionally -- `string(default)`
# to build the fallback, then `strip` to produce a SubString, then a parse -- for
# a measured 432 bytes per call. This function is on the RHS hot path (the flat
# queue's cost model calls it per effector per RHS call, and
# _rhs_effectors_have_heavy_or_heterogeneous_cost calls it per effector inside
# the routing decision, which itself re-runs every RHS call unless the plan step
# cache is on), so at constellation scale that reached multiple GB of pure
# garbage per solve. The value is already snapshotted per run in
# RhsPlanEnvConfig.effector_cost_ns_per_item_default -- this just reads it from
# there. Same hoist as commit d4f5cce4 did for the other per-step env reads;
# these cost-model ones were missed.
@inline function _rhs_effector_static_cost_ns(effector, cost_ns_per_item_default::Float64)::Float64
    rank = _rhs_effector_cost_rank(effector)
    return cost_ns_per_item_default * Float64(rank * rank)
end

@inline function _rhs_effector_estimated_cost_ns(shared_buffers, dynamic_effectors::Tuple, eff_idx::Int)::Float64
    cost_default = _rhs_env_config_from_buffers(shared_buffers).effector_cost_ns_per_item_default
    fallback = _rhs_effector_static_cost_ns(dynamic_effectors[eff_idx], cost_default)
    return _rhs_effector_observed_cost_ns(shared_buffers, eff_idx, fallback)
end

"""
    ConstellationExecutionPlan

Execution-plan metadata for constellation-scale RHS work. The hot path still stores
node work as encoded integers in `SharedBuffers` to avoid per-step allocations; this
type names that layout and leaves a first-class slot for future satellite-satellite
interaction edges.
"""
Base.@kwdef struct ConstellationExecutionPlan
    num_sats::Int
    n_effectors::Int
    workers::Int
    node_count::Int
    edge_count::Int = 0
    partition::Union{Nothing, Symbol} = nothing
    use_packets::Bool = false
end

struct ConstellationInteractionEdgeWorkItem
    source_sat_idx::Int
    target_sat_idx::Int
    eff_idx::Int
end

@inline function _constellation_node_work_item(sat_idx::Int, eff_idx::Int, n_effectors::Int)::Int
    return (sat_idx - 1) * n_effectors + eff_idx
end

@inline function _constellation_node_sat_idx(item::Int, n_effectors::Int)::Int
    return fld(item - 1, n_effectors) + 1
end

@inline function _constellation_node_eff_idx(item::Int, n_effectors::Int)::Int
    return mod(item - 1, n_effectors) + 1
end

# Per-effector "does this belong in the flat queue" mask, built by peeling the
# effector tuple one concretely typed element at a time so every predicate call
# is a static dispatch. Returns NTuple{N, Bool}, which callers can index at a
# runtime eff_idx without allocating (unlike the heterogeneous effector tuple
# itself). Mirrors the skip logic that used to live inline in
# _prepare_rhs_flat_work_items!'s inner loop.
@inline _flat_selection_mask(::Tuple{}, partition)::Tuple{} = ()

@inline function _flat_selection_mask(effs::Tuple, partition)
    effector = first(effs)
    keep = _flat_partition_selected(effector, partition) &&
        # Effectors resolved by pre-passes already wrote into totals.
        !(partition === nothing && (_batchable_effector(effector) || _harmonics_prepass_effector(effector)))
    return (keep, _flat_selection_mask(Base.tail(effs), partition)...)
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

    # Which effectors belong in the flat queue depends only on eff_idx and the
    # partition -- never on sat_idx -- so resolve it once, outside the per-
    # satellite loop, into a homogeneous NTuple{N, Bool}.
    #
    # The previous form indexed `dynamic_effectors[eff_idx]` inside the inner
    # loop, i.e. num_sats * n_effectors times per RHS call. That tuple is
    # heterogeneous, so a runtime index infers to the Union of its element types
    # and allocates a measured 144 bytes each time (0 bytes for a homogeneous
    # tuple), which made this function one of the largest allocation sites in the
    # 12-thread profile of the 1024-satellite atmosphere case. Indexing the Bool
    # mask instead is allocation-free because NTuple{N, Bool} is homogeneous, and
    # building the mask costs n_effectors type-stable predicate calls per RHS
    # call rather than num_sats * n_effectors boxed ones.
    selected = _flat_selection_mask(dynamic_effectors, partition)
    count_items = 0
    @inbounds for sat_idx in 1:num_sats
        p.is_active[sat_idx] || continue
        for eff_idx in 1:n_effectors
            selected[eff_idx] || continue
            count_items += 1
            work_items[count_items] = _constellation_node_work_item(sat_idx, eff_idx, n_effectors)
        end
    end
    # Cost only varies by eff_idx, not sat_idx, and the pre-pass filter above already
    # strips out batchable/harmonics effectors -- so once at most one effector type
    # remains in the flat queue (the common case once a constellation's gravity model
    # is pre-passed away, leaving e.g. just aero), every item shares the same cost key
    # and sorting is a guaranteed no-op that still pays full O(n log n) comparison cost,
    # each comparison re-deriving the same handful of per-effector cost lookups. Skip
    # it entirely in that case; only heterogeneous multi-effector residual queues
    # benefit from (and need) the sort.
    if count_items > 1 && _count_flat_queue_only_effectors(dynamic_effectors) > 1
        # Item cost only depends on eff_idx, so there are just n_effectors distinct
        # values. Resolve each effector's estimated cost once up front; the sort
        # comparator then reads a small tuple instead of re-deriving the cost-model
        # lookup O(n log n) times inside sort!.
        eff_costs = ntuple(
            eff_idx -> _rhs_effector_estimated_cost_ns(p.shared_buffers, dynamic_effectors, eff_idx),
            length(dynamic_effectors),
        )
        sort!(
            @view(work_items[1:count_items]);
            by=item -> -eff_costs[_constellation_node_eff_idx(item, n_effectors)],
            alg=Base.Sort.DEFAULT_STABLE,
        )
    end
    return count_items
end

@inline function _rhs_flat_item_eff_idx(item::Int, n_effectors::Int)::Int
    return _constellation_node_eff_idx(item, n_effectors)
end

@inline function _rhs_flat_item_estimated_cost_ns(shared_buffers, dynamic_effectors::Tuple, item::Int)::Float64
    eff_idx = _rhs_flat_item_eff_idx(item, length(dynamic_effectors))
    return _rhs_effector_estimated_cost_ns(shared_buffers, dynamic_effectors, eff_idx)
end

function _build_constellation_execution_plan!(
    work_items::Vector{Int},
    p,
    dynamic_effectors::Tuple,
    num_sats::Int,
    workers::Int,
    partition::Union{Nothing, Symbol},
)::ConstellationExecutionPlan
    n_effectors = length(dynamic_effectors)
    node_count = _prepare_rhs_flat_work_items!(
        work_items,
        p,
        dynamic_effectors,
        num_sats,
        partition,
    )
    return ConstellationExecutionPlan(
        num_sats=num_sats,
        n_effectors=n_effectors,
        workers=workers,
        node_count=node_count,
        edge_count=0,
        partition=partition,
        use_packets=false,
    )
end

@inline function _with_packet_scheduler(exec_plan::ConstellationExecutionPlan, use_packets::Bool)::ConstellationExecutionPlan
    return ConstellationExecutionPlan(
        num_sats=exec_plan.num_sats,
        n_effectors=exec_plan.n_effectors,
        workers=exec_plan.workers,
        node_count=exec_plan.node_count,
        edge_count=exec_plan.edge_count,
        partition=exec_plan.partition,
        use_packets=use_packets,
    )
end

function _prepare_rhs_flat_work_packets!(
    packet_starts::Vector{Int},
    packet_ends::Vector{Int},
    packet_costs::Vector{Float64},
    packet_elapsed_ns::Vector{Int64},
    work_items::Vector{Int},
    count_items::Int,
    shared_buffers,
    dynamic_effectors::Tuple,
    workers::Int,
)::Int
    count_items <= 0 && return 0
    if length(packet_starts) < count_items
        resize!(packet_starts, count_items)
        resize!(packet_ends, count_items)
        resize!(packet_costs, count_items)
        resize!(packet_elapsed_ns, count_items)
    end

    total_cost_ns = 0.0
    @inbounds for idx in 1:count_items
        total_cost_ns += _rhs_flat_item_estimated_cost_ns(shared_buffers, dynamic_effectors, work_items[idx])
    end
    target_ns = max(_rhs_env_config_from_buffers(shared_buffers).flat_packet_target_min_ns, total_cost_ns / max(1, workers * 4))

    packet_count = 0
    idx = 1
    @inbounds while idx <= count_items
        packet_count += 1
        packet_starts[packet_count] = idx
        packet_elapsed_ns[packet_count] = Int64(0)
        running_cost = 0.0
        first_item = true
        while idx <= count_items
            item_cost = _rhs_flat_item_estimated_cost_ns(shared_buffers, dynamic_effectors, work_items[idx])
            if !first_item && running_cost + item_cost > target_ns
                break
            end
            running_cost += item_cost
            idx += 1
            first_item = false
            if running_cost >= target_ns
                break
            end
        end
        packet_ends[packet_count] = idx - 1
        packet_costs[packet_count] = max(running_cost, 1.0)
    end
    return packet_count
end

function _rhs_flat_packet_work_stats(
    shared_buffers,
    dynamic_effectors::Tuple,
    work_items::Vector{Int},
    count_items::Int,
)::Tuple{Float64, Float64}
    total_cost_ns = 0.0
    min_cost_ns = Inf
    max_cost_ns = 0.0
    @inbounds for idx in 1:count_items
        item_cost = _rhs_flat_item_estimated_cost_ns(shared_buffers, dynamic_effectors, work_items[idx])
        total_cost_ns += item_cost
        min_cost_ns = min(min_cost_ns, item_cost)
        max_cost_ns = max(max_cost_ns, item_cost)
    end
    heterogeneity = max_cost_ns / max(min_cost_ns, 1.0)
    return total_cost_ns, heterogeneity
end

function _rhs_flat_use_packet_scheduler(
    shared_buffers,
    dynamic_effectors::Tuple,
    work_items::Vector{Int},
    count_items::Int,
)::Bool
    env = _rhs_env_config_from_buffers(shared_buffers)
    mode = env.flat_packet_scheduler_mode
    mode == :off && return false
    if mode == :on
        return count_items > 1
    end
    if hasproperty(shared_buffers, :rhs_flat_packet_disabled) && shared_buffers.rhs_flat_packet_disabled[]
        return false
    end
    count_items >= env.flat_packet_min_items || return false
    total_cost_ns, heterogeneity = _rhs_flat_packet_work_stats(shared_buffers, dynamic_effectors, work_items, count_items)
    total_cost_ns >= env.flat_packet_work_ns_threshold || return false
    heterogeneity >= env.flat_packet_heterogeneity_threshold || return false
    return true
end

function _update_rhs_flat_packet_overhead_model!(
    shared_buffers,
    packet_overhead_ns::Int64,
    flat_elapsed_ns::Int64,
)::Nothing
    packet_overhead_ns > 0 || return nothing
    flat_elapsed_ns > 0 || return nothing
    if !(hasproperty(shared_buffers, :rhs_flat_packet_overhead_ema) &&
         hasproperty(shared_buffers, :rhs_flat_packet_overhead_samples) &&
         hasproperty(shared_buffers, :rhs_flat_packet_disabled))
        return nothing
    end
    env = _rhs_env_config_from_buffers(shared_buffers)
    ratio = min(1.0, Float64(packet_overhead_ns) / Float64(flat_elapsed_ns))
    α = env.effector_cost_ema_alpha
    previous = shared_buffers.rhs_flat_packet_overhead_ema[]
    shared_buffers.rhs_flat_packet_overhead_ema[] = if isfinite(previous) && previous >= 0.0
        (1.0 - α) * previous + α * ratio
    else
        ratio
    end
    shared_buffers.rhs_flat_packet_overhead_samples[] =
        min(typemax(Int64), shared_buffers.rhs_flat_packet_overhead_samples[] + Int64(1))
    if shared_buffers.rhs_flat_packet_overhead_samples[] >= env.flat_packet_overhead_min_samples &&
       shared_buffers.rhs_flat_packet_overhead_ema[] >= env.flat_packet_overhead_disable_ratio
        shared_buffers.rhs_flat_packet_disabled[] = true
    end
    return nothing
end

function _update_rhs_flat_packet_cost_model!(
    shared_buffers,
    dynamic_effectors::Tuple,
    work_items::Vector{Int},
    packet_starts::Vector{Int},
    packet_ends::Vector{Int},
    packet_costs::Vector{Float64},
    packet_elapsed_ns::Vector{Int64},
    packet_count::Int,
)::Nothing
    packet_count <= 0 && return nothing
    n_effectors = length(dynamic_effectors)
    @inbounds for packet_idx in 1:packet_count
        elapsed_ns = packet_elapsed_ns[packet_idx]
        elapsed_ns > 0 || continue
        estimated_packet_cost = max(packet_costs[packet_idx], 1.0)
        scale = Float64(elapsed_ns) / estimated_packet_cost
        for item_idx in packet_starts[packet_idx]:packet_ends[packet_idx]
            item = work_items[item_idx]
            eff_idx = _rhs_flat_item_eff_idx(item, n_effectors)
            item_estimate = _rhs_effector_estimated_cost_ns(shared_buffers, dynamic_effectors, eff_idx)
            _update_rhs_effector_cost_model!(shared_buffers, eff_idx, max(1.0, item_estimate * scale))
        end
    end
    return nothing
end

# ── Batchable-effector trait ────────────────────────────────────────────────────
# Effectors that can be evaluated over all satellites in a single pre-pass,
# reading shared-environment data (ephemeris, sun position) that was prefilled
# once by _prefill_shared_body_samples!.  The pre-pass writes directly into the
# totals matrix, and those effectors are skipped in the per-(sat,effector) flat
# queue.  This eliminates per-item channel/worker dispatch overhead for NBody,
# SRP, and inverse-square gravity, which is the dominant cost for high-thread-count
# constellations once the ephemeris cache is warm (NBody/SRP) or the whole point
# for a kernel too cheap to ever amortise per-item dispatch (inverse-square gravity).
#
# GRAM/aero are NOT batchable: density is altitude/location-dependent per
# satellite and mutable model state prevents safe cross-satellite sharing.
#
# Inverse-square gravity is only batchable when it produces no torque: the batch
# kernel writes force only, so a model configured with `gravity_gradient=true`
# (which needs the per-satellite quaternion) must keep going through the
# per-satellite wrench path instead.
@inline _batchable_effector(::Any)::Bool = false
@inline _batchable_effector(::SimulationModel.NBodyGravityModel)::Bool = true
@inline _batchable_effector(::SimulationModel.SolarRadiationPressureModel)::Bool = true
@inline _batchable_effector(effector::SimulationModel.InverseSquaredGravityModel)::Bool = !effector.gravity_gradient
@inline _batchable_effector(effector::SimulationModel.InverseSquaredJ2GravityModel)::Bool = !effector.gravity_gradient

@inline function _has_any_batchable_effector(effectors::Tuple)::Bool
    @inbounds for e in effectors
        _batchable_effector(e) && return true
    end
    return false
end

@inline function _count_non_batchable_effectors(effectors::Tuple)::Int
    count = 0
    @inbounds for e in effectors
        _batchable_effector(e) || (count += 1)
    end
    return count
end

# Harmonics has its own SIMD pre-pass (different from the batchable trait because
# it uses mutable LPI matrix state and requires sc_state rather than pos_buffers).
@inline _harmonics_prepass_effector(::Any)::Bool = false
@inline _harmonics_prepass_effector(::SimulationModel.GravitationalHarmonicsModel)::Bool = true

@inline function _has_any_harmonics_effector(effectors::Tuple)::Bool
    @inbounds for e in effectors
        _harmonics_prepass_effector(e) && return true
    end
    return false
end

# Effectors that still need flat queue dispatch (not handled by any pre-pass).
@inline function _count_flat_queue_only_effectors(effectors::Tuple)::Int
    count = 0
    @inbounds for e in effectors
        _batchable_effector(e) || _harmonics_prepass_effector(e) || (count += 1)
    end
    return count
end

# NBody batch pre-pass: third-body positions are read from the SpiceRhsMemo (or
# the pre-warmed NBodyEphemerisCache) exactly once, then the force is accumulated
# for all active satellites with a tight scalar loop.  No allocations.
function _accumulate_nbody_flat_batch!(
    totals::Matrix{Float64},
    effector::SimulationModel.NBodyGravityModel,
    pos_buffers::Vector{SVector{3, Float64}},
    mass_buffers::Vector{Float64},
    active_flags,
    p,
    t::Float64,
    num_sats::Int,
)::Nothing
    n_bodies = length(effector.body_names)
    n_bodies == 0 && return nothing

    et = p.shared_buffers.et_start[] + t
    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(effector.primary_body_name)
    memo_enabled  = p.shared_buffers.spice_rhs_memo_enabled[]
    memo          = p.shared_buffers.spice_rhs_memo
    nbody_cache   = p.shared_buffers.nbody_ephemeris_cache[]
    counter       = p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls
    pert          = SimulationModel.DynamicEffectors.PerturbationEffectors

    # Gather third-body positions once (O(n_bodies), not O(n_bodies × n_sats)).
    body_positions = ntuple(n_bodies) do k
        bname = SimulationModel.DynamicEffectors._spice_query_name(effector.body_names[k])
        if nbody_cache isa SimulationModel.NBodyEphemerisCache
            pos = SimulationModel.DynamicEffectors._nbody_body_position_from_cache_j2000_m(
                nbody_cache, et, bname, primary_body_name)
            pos !== nothing && return pos
        end
        pert._nbody_body_position_from_spice_j2000_m(bname, et, primary_body_name, memo_enabled, memo, counter)
    end

    body_mus = effector.body_mus
    @inbounds for sat_idx in 1:num_sats
        active_flags[sat_idx] || continue
        r = pos_buffers[sat_idx]
        mass = mass_buffers[sat_idx]
        r1, r2, r3 = r[1], r[2], r[3]
        fx, fy, fz = 0.0, 0.0, 0.0
        for k in 1:n_bodies
            rk = body_positions[k]
            rk1, rk2, rk3 = rk[1], rk[2], rk[3]
            dx = rk1 - r1; dy = rk2 - r2; dz = rk3 - r3
            d2 = dx*dx + dy*dy + dz*dz
            d3 = d2 * sqrt(d2)
            rk2_norm = rk1*rk1 + rk2*rk2 + rk3*rk3
            rk3_norm = rk2_norm * sqrt(rk2_norm)
            mu = body_mus[k]
            fx += mu * (dx/d3 - rk1/rk3_norm) * mass
            fy += mu * (dy/d3 - rk2/rk3_norm) * mass
            fz += mu * (dz/d3 - rk3/rk3_norm) * mass
        end
        totals[1, sat_idx] += fx
        totals[2, sat_idx] += fy
        totals[3, sat_idx] += fz
        # torques[4..6] stay zero (NBody produces no torque)
    end
    return nothing
end

# SRP batch pre-pass: sun position is read once (from ephemeris cache or memo),
# then SRP acceleration is accumulated for all active satellites.
function _accumulate_srp_flat_batch!(
    totals::Matrix{Float64},
    effector::SimulationModel.SolarRadiationPressureModel,
    pos_buffers::Vector{SVector{3, Float64}},
    mass_buffers::Vector{Float64},
    active_flags,
    p,
    t::Float64,
    num_sats::Int,
)::Nothing
    effector.A == 0.0 && return nothing

    et = p.shared_buffers.et_start[] + t
    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(p.args.environment_model.planet.name)
    memo_enabled = p.shared_buffers.spice_rhs_memo_enabled[]
    memo         = p.shared_buffers.spice_rhs_memo
    counter      = p.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls
    pert         = SimulationModel.DynamicEffectors.PerturbationEffectors

    sun_pos = if effector.direct || effector.albedo
        srp_cache = p.shared_buffers.srp_sun_ephemeris_cache[]
        if srp_cache isa SimulationModel.SRPSunEphemerisCache
            pos = SimulationModel.DynamicEffectors._srp_sun_position_from_cache_j2000_m(srp_cache, et)
            pos !== nothing ? pos : pert._srp_sun_position_from_spice_j2000_m(et, primary_body_name, memo_enabled, memo, counter)
        else
            pert._srp_sun_position_from_spice_j2000_m(et, primary_body_name, memo_enabled, memo, counter)
        end
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    planet = p.args.environment_model.planet
    @inbounds for sat_idx in 1:num_sats
        active_flags[sat_idx] || continue
        pos_ii = pos_buffers[sat_idx]
        mass   = mass_buffers[sat_idx]
        accel  = pert._srp_total_acceleration_ii(effector, planet, pos_ii, sun_pos, mass)
        totals[1, sat_idx] += mass * accel[1]
        totals[2, sat_idx] += mass * accel[2]
        totals[3, sat_idx] += mass * accel[3]
        # torques[4..6] stay zero
    end
    return nothing
end

# Inverse-square gravity batch pre-pass: no shared ephemeris lookup to hoist (the
# primary is at the frame origin), so the win here is purely eliminating
# per-satellite Polyester/flat-queue dispatch for a kernel whose body is a few
# FLOPs — too cheap to ever amortise task-spawn overhead. Delegates to the same
# `_inverse_squared_gravity_accel` helper used by `calcForceTorque`/`wrench`, so
# results are bit-identical to the per-satellite path.
function _accumulate_invsq_flat_batch!(
    totals::Matrix{Float64},
    effector::SimulationModel.InverseSquaredGravityModel,
    pos_buffers::Vector{SVector{3, Float64}},
    mass_buffers::Vector{Float64},
    active_flags,
    p,
    t::Float64,
    num_sats::Int,
)::Nothing
    planet = p.args.environment_model.planet
    grav = SimulationModel.DynamicEffectors.GravityEffectors
    @inbounds for sat_idx in 1:num_sats
        active_flags[sat_idx] || continue
        pos_ii = pos_buffers[sat_idx]
        mass   = mass_buffers[sat_idx]
        accel  = grav._inverse_squared_gravity_accel(pos_ii, planet)
        totals[1, sat_idx] += mass * accel[1]
        totals[2, sat_idx] += mass * accel[2]
        totals[3, sat_idx] += mass * accel[3]
        # torques[4..6] stay zero (only reached when !effector.gravity_gradient)
    end
    return nothing
end

# Inverse-square + J2 gravity batch pre-pass: J2 needs the position in the
# planet-fixed frame, so l_pi (inertial->planet rotation) is computed once per
# RHS call — analogous to the harmonics pre-pass computing lpi once instead of
# per satellite — then the per-satellite loop is a tight, allocation-free scalar
# sweep with no further shared lookups.
function _accumulate_invsq_j2_flat_batch!(
    totals::Matrix{Float64},
    effector::SimulationModel.InverseSquaredJ2GravityModel,
    pos_buffers::Vector{SVector{3, Float64}},
    mass_buffers::Vector{Float64},
    active_flags,
    p,
    t::Float64,
    num_sats::Int,
)::Nothing
    planet = p.args.environment_model.planet
    l_pi = _planet_lpi_at_engine(p, t)
    grav = SimulationModel.DynamicEffectors.GravityEffectors
    @inbounds for sat_idx in 1:num_sats
        active_flags[sat_idx] || continue
        pos_ii = pos_buffers[sat_idx]
        mass   = mass_buffers[sat_idx]
        pos_pp = SVector{3, Float64}(l_pi * pos_ii)
        accel_pp = grav._inverse_squared_j2_gravity_accel(pos_pp, planet)
        accel_ii = l_pi' * accel_pp
        totals[1, sat_idx] += mass * accel_ii[1]
        totals[2, sat_idx] += mass * accel_ii[2]
        totals[3, sat_idx] += mass * accel_ii[3]
        # torques[4..6] stay zero (only reached when !effector.gravity_gradient)
    end
    return nothing
end

# Dispatch to the appropriate batch kernel for a single batchable effector.
@inline function _accumulate_batchable_effector_flat!(
    totals::Matrix{Float64},
    effector,
    pos_buffers::Vector{SVector{3, Float64}},
    mass_buffers::Vector{Float64},
    active_flags,
    p,
    t::Float64,
    num_sats::Int,
)::Nothing
    if effector isa SimulationModel.NBodyGravityModel
        return _accumulate_nbody_flat_batch!(totals, effector, pos_buffers, mass_buffers, active_flags, p, t, num_sats)
    elseif effector isa SimulationModel.SolarRadiationPressureModel
        return _accumulate_srp_flat_batch!(totals, effector, pos_buffers, mass_buffers, active_flags, p, t, num_sats)
    elseif effector isa SimulationModel.InverseSquaredGravityModel
        return _accumulate_invsq_flat_batch!(totals, effector, pos_buffers, mass_buffers, active_flags, p, t, num_sats)
    elseif effector isa SimulationModel.InverseSquaredJ2GravityModel
        return _accumulate_invsq_j2_flat_batch!(totals, effector, pos_buffers, mass_buffers, active_flags, p, t, num_sats)
    end
    return nothing
end

# Forward a plan's `scheduler` to the dispatch primitives only when the plan came
# from pre-solve calibration.
#
# Heuristic plans in setup.jl have carried a `scheduler` field since long before
# the dispatch primitives honoured one -- they always read
# SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER instead -- so those values have
# never been measured against anything. Honouring them now would silently change
# R0-R3 too, which is a different change from making the adaptive profiles'
# scheduler a routed decision. A calibrated plan's scheduler is a swept, timed
# choice, so that one is forwarded; everything else resolves to `:auto`, i.e.
# the env var, exactly as before.
#
# The override Ref is also what the calibration sweep writes each candidate
# into, so the sweep measures each candidate under its own scheduler.
@inline function _dispatch_scheduler(p, plan)::Symbol
    if p !== nothing && hasproperty(p, :shared_buffers) &&
       p.shared_buffers.rhs_plan_override[] !== nothing
        return plan.scheduler
    end
    return :auto
end

function _accumulate_harmonics_flat_batch!(
    sc_state,
    p,
    t::Float64,
    model::SimulationModel.GravitationalHarmonicsModel,
    plan;
    init_scratch::Bool=true,
)::Nothing
    num_sats = length(sc_state)
    active_sats = count(identity, p.is_active)
    active_sats <= 0 && return nothing
    rhs_env = _rhs_env_config(p)
    min_sats = rhs_env.harmonics_batch_spin_barrier ?
        1 : rhs_env.harmonics_batch_min_sats_per_worker
    capped_allotment = max(1, min(plan.allotment, fld(active_sats, max(1, min_sats))))
    workers = SimulationModel.ParallelPolicy.thread_worker_count(active_sats, capped_allotment)
    if init_scratch
        _ensure_rhs_flat_effector_scratch!(p.shared_buffers, num_sats, workers; zero_partials=false)
    end
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
    pool = SimulationModel.DynamicEffectors.PerturbationEffectors._get_harmonics_batch_pool_cached!(
        p.shared_buffers.rhs_harmonics_batch_pool, model, n_workers, batch_size,
    )
    if n_workers <= 1
        SimulationModel.DynamicEffectors.PerturbationEffectors._harmonics_flat_batch_kernel!(
            totals, model, sc_state, work_items, 1, count_items, lpi, pool[1]
        )
    else
        dispatch_fn = rhs_env.harmonics_batch_spin_barrier ?
            SimulationModel.ParallelPolicy.threaded_foreach_worker_spin :
            SimulationModel.ParallelPolicy.threaded_foreach_worker_persistent
        dispatch_fn(
            :rhs_harmonics_batch,
            n_workers,
            plan.allotment;
            scheduler=_dispatch_scheduler(p, plan),
        ) do _worker_id, w
            item_start = (w - 1) * batch_size + 1
            item_end   = min(w * batch_size, count_items)
            item_start > count_items && return
            SimulationModel.DynamicEffectors.PerturbationEffectors._harmonics_flat_batch_kernel!(
                totals, model, sc_state, work_items, item_start, item_end, lpi, pool[w]
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
            env=_policy_env_config(p),
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
    # Partials are only written by the flat queue itself; pre-pass-only effector
    # sets (e.g. harmonics-only constellations, which take the batch-kernel
    # shortcut below and write straight into totals) never touch them, so skip
    # the O(6·N·W) zeroing sweep in that case. Partitioned calls always run the
    # queue for their selected effectors.
    needs_flat_queue = partition === nothing ?
        _count_flat_queue_only_effectors(dynamic_effectors) > 0 : true
    _ensure_rhs_flat_effector_scratch!(p.shared_buffers, num_sats, workers; zero_partials=needs_flat_queue)
    selected_count = partition === nothing ? n_effectors : _partition_selected_count(dynamic_effectors, partition)
    selected_count <= 0 && return nothing

    if partition === nothing &&
       n_effectors == 1 &&
       dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel
        return _accumulate_harmonics_flat_batch!(sc_state, p, t, dynamic_effectors[1], plan)
    end

    partials = p.shared_buffers.rhs_flat_effector_partials[]
    totals = p.shared_buffers.rhs_flat_effector_totals[]
    work_items = p.shared_buffers.rhs_flat_work_items[]
    packet_starts = p.shared_buffers.rhs_flat_packet_starts[]
    packet_ends = p.shared_buffers.rhs_flat_packet_ends[]
    packet_costs = p.shared_buffers.rhs_flat_packet_costs[]
    packet_elapsed_ns = p.shared_buffers.rhs_flat_packet_elapsed_ns[]
    _ensure_rhs_effector_cost_model!(p.shared_buffers, n_effectors)

    needs_state_sample = if partition === nothing
        any(_wrench_method_available(effector) for effector in dynamic_effectors) ||
            _has_any_batchable_effector(dynamic_effectors)
    else
        _partition_needs_state_sample(dynamic_effectors, partition)
    end
    if needs_state_sample
        _prefill_rhs_flat_state_samples!(p.shared_buffers, sc_state, p)
    end
    spacecraft = p.args.dynamics_model.spacecraft
    orientation_sim = p.args.mission_configuration.orientation_sim

    # ── Batchable effector pre-pass (NBody, SRP) ─────────────────────────────
    # Run batchable effectors as a single serial sweep before the flat queue.
    # They read shared environment data (third-body/sun positions) that was
    # prefilled once by _prefill_shared_body_samples!, avoiding per-item channel
    # dispatch overhead.  Only applies to the unpartitioned path; IMEX-partitioned
    # calls let the flat queue handle partition filtering.
    if partition === nothing && _has_any_batchable_effector(dynamic_effectors)
        pos_buffers  = p.shared_buffers.rhs_flat_state_pos_ii[]
        mass_buffers = p.shared_buffers.rhs_flat_state_mass_kg[]
        @inbounds for effector in dynamic_effectors
            _batchable_effector(effector) || continue
            _accumulate_batchable_effector_flat!(totals, effector, pos_buffers, mass_buffers, p.is_active, p, t, num_sats)
        end
        _count_non_batchable_effectors(dynamic_effectors) == 0 && return nothing
    end

    # ── Harmonics SIMD pre-pass ───────────────────────────────────────────────
    # Run the SIMD batch kernel for GravitationalHarmonicsModel as a pre-pass,
    # writing directly into the already-zeroed totals matrix.  Skips scratch
    # re-initialization so batchable contributions above are preserved.  The
    # single-effector shortcut at the top of this function already handles the
    # harmonics-only case; this pre-pass covers the multi-effector path.
    if partition === nothing && _has_any_harmonics_effector(dynamic_effectors)
        @inbounds for effector in dynamic_effectors
            effector isa SimulationModel.GravitationalHarmonicsModel || continue
            _accumulate_harmonics_flat_batch!(sc_state, p, t, effector, plan; init_scratch=false)
            break
        end
        _count_flat_queue_only_effectors(dynamic_effectors) == 0 && return nothing
    end

    needs_timing = plan.policy_applied
    started_ns = needs_timing ? time_ns() : UInt64(0)
    exec_plan = _build_constellation_execution_plan!(
        work_items,
        p,
        dynamic_effectors,
        num_sats,
        workers,
        partition,
    )
    count_items = exec_plan.node_count
    count_items <= 0 && return nothing

    use_packets = _rhs_flat_use_packet_scheduler(p.shared_buffers, dynamic_effectors, work_items, count_items)
    exec_plan = _with_packet_scheduler(exec_plan, use_packets)
    packet_count = 0
    packet_overhead_ns = Int64(0)
    if exec_plan.use_packets
        packet_prepare_started_ns = needs_timing ? time_ns() : UInt64(0)
        packet_count = _prepare_rhs_flat_work_packets!(
            packet_starts,
            packet_ends,
            packet_costs,
            packet_elapsed_ns,
            work_items,
            count_items,
            p.shared_buffers,
            dynamic_effectors,
            exec_plan.workers,
        )
        packet_count <= 0 && return nothing
        if needs_timing
            packet_overhead_ns += Int64(time_ns() - packet_prepare_started_ns)
        end

        SimulationModel.ParallelPolicy.threaded_foreach_worker_persistent(
            :rhs_flat_queue_packets, packet_count, plan.allotment; scheduler=_dispatch_scheduler(p, plan)
        ) do worker_id, packet_idx
            packet_started_ns = needs_timing ? time_ns() : UInt64(0)
            @inbounds for item_idx in packet_starts[packet_idx]:packet_ends[packet_idx]
                item = work_items[item_idx]
                sat_idx = _constellation_node_sat_idx(item, exec_plan.n_effectors)
                eff_idx = _constellation_node_eff_idx(item, exec_plan.n_effectors)
                effector = dynamic_effectors[eff_idx]
                partition === nothing && (_batchable_effector(effector) || _harmonics_prepass_effector(effector)) && continue
                @views sc_view = sc_state[sat_idx]
                state_sample = _wrench_method_available(effector) ?
                    _rhs_flat_state_sample_from_buffers(p.shared_buffers, spacecraft, sat_idx, orientation_sim) :
                    nothing
                force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
                partials[1, sat_idx, worker_id] += force[1]
                partials[2, sat_idx, worker_id] += force[2]
                partials[3, sat_idx, worker_id] += force[3]
                partials[4, sat_idx, worker_id] += torque[1]
                partials[5, sat_idx, worker_id] += torque[2]
                partials[6, sat_idx, worker_id] += torque[3]
            end
            if needs_timing
                @inbounds packet_elapsed_ns[packet_idx] = Int64(time_ns() - packet_started_ns)
            end
            return nothing
        end
    else
        SimulationModel.ParallelPolicy.threaded_foreach_worker_persistent(
            :rhs_flat_queue, count_items, plan.allotment; scheduler=_dispatch_scheduler(p, plan)
        ) do worker_id, item_idx
            item = work_items[item_idx]
            sat_idx = _constellation_node_sat_idx(item, exec_plan.n_effectors)
            eff_idx = _constellation_node_eff_idx(item, exec_plan.n_effectors)
            effector = dynamic_effectors[eff_idx]
            partition === nothing && (_batchable_effector(effector) || _harmonics_prepass_effector(effector)) && return nothing
            @views sc_view = sc_state[sat_idx]
            state_sample = _wrench_method_available(effector) ?
                _rhs_flat_state_sample_from_buffers(p.shared_buffers, spacecraft, sat_idx, orientation_sim) :
                nothing
            force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
            partials[1, sat_idx, worker_id] += force[1]
            partials[2, sat_idx, worker_id] += force[2]
            partials[3, sat_idx, worker_id] += force[3]
            partials[4, sat_idx, worker_id] += torque[1]
            partials[5, sat_idx, worker_id] += torque[2]
            partials[6, sat_idx, worker_id] += torque[3]
            return nothing
        end
    end

    # Worker-major reduction: partials is column-major 6×N×W, so for a fixed
    # worker the satellite dimension is contiguous (stride 6); satellite-major
    # order with the worker innermost would jump 6N doubles per iteration.
    @inbounds for worker_id in 1:exec_plan.workers
        for sat_idx in 1:num_sats
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
            env=_policy_env_config(p),
        )
        if exec_plan.use_packets
            feedback_started_ns = time_ns()
            _update_rhs_flat_packet_cost_model!(
                p.shared_buffers,
                dynamic_effectors,
                work_items,
                packet_starts,
                packet_ends,
                packet_costs,
                packet_elapsed_ns,
                packet_count,
            )
            packet_overhead_ns += Int64(time_ns() - feedback_started_ns)
            _update_rhs_flat_packet_overhead_model!(p.shared_buffers, packet_overhead_ns, elapsed_ns)
        end
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
    p.shared_buffers.rhs_solar_prefilled[] = false

    for effector in dynamic_effectors
        needs_solar = effector isa SimulationModel.SolarRadiationPressureModel ||
            (_wrench_method_available(effector) && SimulationModel.environment_requirements(effector).solar)
        if needs_solar
            solar = sample_solar_ephemeris(sc_view, p, first_active, t)
            p.shared_buffers.rhs_flat_solar_pos_ii[] = solar.sun_pos_ii
            p.shared_buffers.rhs_flat_solar_t[] = t
            p.shared_buffers.rhs_solar_prefilled[] = true
            break
        end
    end

    for effector in dynamic_effectors
        if effector isa SimulationModel.NBodyGravityModel
            sample_third_body_ephemerides(effector, sc_view, p, first_active, t)
            break
        end
    end

    for effector in dynamic_effectors
        if effector isa SimulationModel.GravitationalHarmonicsModel
            et = p.shared_buffers.et_start[] + t
            SimulationModel.DynamicEffectors.PerturbationEffectors._harmonics_lpi_at!(effector, p, et)
            break
        end
    end
    return nothing
end

# Parallel environment pre-sample: fills planet-frame component arrays for all active
# satellites and, when requested, shared_buffers.densities/temperatures/winds before the
# flat effector queue. l_pi is pre-computed once so the batch avoids acquiring
# harmonics_lpi_lock per satellite.
function _prefill_environment_samples!(p, t::Float64, sc_state; atmosphere::Bool)::Nothing
    num_sats = length(sc_state)
    num_sats == 0 && return nothing
    planet = p.args.environment_model.planet
    l_pi = _planet_lpi_at_engine(p, t)
    p.shared_buffers.rhs_flat_planet_lpi[] = l_pi
    planet_pos = p.shared_buffers.rhs_flat_planet_pos_pp[]
    planet_vel = p.shared_buffers.rhs_flat_planet_vel_pp[]
    planet_alt = p.shared_buffers.rhs_flat_planet_alt_m[]
    planet_lat = p.shared_buffers.rhs_flat_planet_lat_rad[]
    planet_lon = p.shared_buffers.rhs_flat_planet_lon_rad[]
    if length(planet_pos) < num_sats
        resize!(planet_pos, num_sats)
        resize!(planet_vel, num_sats)
        resize!(planet_alt, num_sats)
        resize!(planet_lat, num_sats)
        resize!(planet_lon, num_sats)
    end

    decision = SimulationModel.SimulationCallbacks._density_callback_thread_decision(p.args, num_sats)
    worker_allotment = decision.use_threads ? decision.allotment : 1

    # Whether the atmosphere half of this pass goes to the distributed density
    # service. Decided once, before the loop, so every satellite in one
    # derivative evaluation takes the same path -- a per-satellite decision would
    # split one batch across two mechanisms and defeat the batching.
    #
    # This is the batch point the service needs, and the reason it can be exact:
    # the loop below already computes every satellite's planet frame before any
    # force is accumulated, so all N density queries for this evaluation are
    # available together, at the true stage state. Nothing is frozen or reused.
    batch_atmosphere = atmosphere &&
        SimulationModel.SimulationCallbacks._rhs_density_service_candidate(p, num_sats)

    SimulationModel.ParallelPolicy.threaded_foreach_worker_persistent(:rhs_atmosphere, num_sats, worker_allotment) do _, sat_idx
        if p.is_active[sat_idx]
            @views sc_view = sc_state[sat_idx]
            planet_frame = sample_planet_frame_with_lpi(sc_view, planet, l_pi)
            @inbounds begin
                planet_pos[sat_idx] = planet_frame.pos_pp
                planet_vel[sat_idx] = planet_frame.vel_pp
                planet_alt[sat_idx] = planet_frame.alt_m
                planet_lat[sat_idx] = planet_frame.lat_rad
                planet_lon[sat_idx] = planet_frame.lon_rad
            end
            if atmosphere && !batch_atmosphere
                _sample_atmosphere_from_planet_frame(sc_view, planet_frame, p, sat_idx, t; write_buffers=true)
            end
        end
    end

    if batch_atmosphere
        served = SimulationModel.SimulationCallbacks._rhs_density_service_fill!(
            p, t, num_sats, planet_alt, planet_lat, planet_lon
        )
        if !served
            # Service declined or a worker failed. The planet frames above are
            # already in the buffers, so the fallback re-reads them rather than
            # recomputing, and the result is identical to never having tried.
            SimulationModel.ParallelPolicy.threaded_foreach_worker_persistent(:rhs_atmosphere, num_sats, worker_allotment) do _, sat_idx
                if p.is_active[sat_idx]
                    @views sc_view = sc_state[sat_idx]
                    _sample_atmosphere_from_planet_frame(
                        sc_view, sample_buffered_planet_frame(p, sat_idx), p, sat_idx, t;
                        write_buffers=true,
                    )
                end
            end
        end
    end
    return nothing
end

function _prefill_atmosphere_samples!(p, t::Float64, sc_state)::Nothing
    return _prefill_environment_samples!(p, t, sc_state; atmosphere=true)
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

    # Phase 2 (parallel): pre-sample reusable per-satellite environment components so
    # wrench effectors in the hot loop read typed buffers instead of rebuilding samples.
    needs_planet_frame_prefill = partition === nothing && any(dynamic_effectors) do eff
        _wrench_method_available(eff) && begin
            req = SimulationModel.environment_requirements(eff)
            req.planet_frame || req.atmosphere
        end
    end
    needs_atm_prefill = partition === nothing && any(dynamic_effectors) do eff
        _wrench_method_available(eff) && SimulationModel.environment_requirements(eff).atmosphere
    end
    if needs_planet_frame_prefill
        _prefill_environment_samples!(p, t, sc_state; atmosphere=needs_atm_prefill)
        p.shared_buffers.rhs_planet_frame_prefilled[] = true
        p.shared_buffers.rhs_atmosphere_prefilled[] = needs_atm_prefill
    end

    try
        _accumulate_dynamic_effectors_flat_batch!(sc_state, p, t, dynamic_effectors, plan; partition=partition)
    finally
        if needs_planet_frame_prefill
            p.shared_buffers.rhs_planet_frame_prefilled[] = false
            p.shared_buffers.rhs_atmosphere_prefilled[] = false
        end
    end
    totals = p.shared_buffers.rhs_flat_effector_totals[]

    SimulationModel.ParallelPolicy.threaded_foreach(length(sc_state), plan.allotment) do i
        @inbounds if !p.is_active[i] ||
                     (rhs_kind == :implicit && _spacecraft_outside_atmosphere_for_current_state(sc_state[i], p, i, t))
            sc_du[i] .= 0.0
            return
        end
        @inbounds @views begin
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                    )
                end
                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            else
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = if rhs_kind == :explicit || rhs_kind == :full
                    _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                else
                    0.0
                end
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
            return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
        elseif hasproperty(sc_state, :vel)
            return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
        end
        return SVector{3, Float64}(sc_state)
    end

    if sat_idx <= length(state)
        sat_state = state[sat_idx]
        if hasproperty(sat_state, :pos)
            return SVector{3, Float64}(sat_state[1], sat_state[2], sat_state[3])
        elseif hasproperty(sat_state, :vel)
            return SVector{3, Float64}(sat_state[1], sat_state[2], sat_state[3])
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
    use_rhs_batch = _rhs_batch_parallel_enabled(p, length(vel_state.sc))

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
    use_rhs_batch = _rhs_batch_parallel_enabled(p, length(q_state))
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
    rw_assembly=nothing,
    rw_torque_body::SVector{3, Float64}=SVector{3, Float64}(0.0, 0.0, 0.0),
)
    omega_body = SimulationModel.DynamicsRotational.body_angular_velocity(sc_view.ω)
    tau_body = SimulationModel.DynamicsRotational.body_torque(torques)

    if propagate_quaternion
        du_view.q .= SimulationModel.DynamicsRotational.quaternion_derivative(omega_body, sc_view.q)
    else
        du_view.q .= 0.0
    end

    h_wheel_body = SVector{3, Float64}(0.0, 0.0, 0.0)
    if rw_assembly !== nothing && rw_assembly.n_wheels > 0
        h_wheel_body = rw_assembly.J_rw * SVector{rw_assembly.n_wheels, Float64}(sc_view.h_wheels)
        du_view.h_wheels .= rw_assembly.J_rw_pinv * (-rw_torque_body)
    end

    du_view.ω .= SimulationModel.DynamicsRotational.angular_acceleration(
        omega_body,
        inertia_tensor,
        tau_body;
        include_gyroscopic=include_gyroscopic,
        h_wheel_body=h_wheel_body,
    )
    return nothing
end

"""Return whether an effector owns a robot-arm plan for the requested spacecraft."""
@inline function _robot_arm_effector_matches(effector, sat_idx::Int)::Bool
    return hasproperty(effector, :plan) &&
        hasproperty(effector, :spacecraft_idx) &&
        getproperty(effector, :spacecraft_idx) == sat_idx &&
        getproperty(effector, :plan) isa SimulationModel.RobotArmPlan
end

"""Read a numeric effector field as `Float64`, falling back to a default value."""
@inline _effector_float(effector, name::Symbol, default::Float64)::Float64 =
    hasproperty(effector, name) ? Float64(getproperty(effector, name)) : default

"""Read an optional effector field, falling back to a default value."""
@inline _effector_value(effector, name::Symbol, default) =
    hasproperty(effector, name) ? getproperty(effector, name) : default

"""Build the cloth-coupled robot-arm RHS configuration from an active effector."""
@inline function _robot_arm_coupling_from_effector(effector, t::Float64)
    return (
        plan=getproperty(effector, :plan),
        t_s=max(0.0, t - _effector_float(effector, :updated_at_s, 0.0)),
        k_translation_n_m=_effector_value(effector, :k_translation_n_m, 5.0e3),
        c_translation_n_s_m=_effector_value(effector, :c_translation_n_s_m, 30.0),
        k_rotation_n_m_rad=_effector_value(effector, :k_rotation_n_m_rad, 15.0),
        c_rotation_n_m_s_rad=_effector_value(effector, :c_rotation_n_m_s_rad, 0.5),
        joint_actuators=hasproperty(effector, :joint_actuators) ?
            getproperty(effector, :joint_actuators) :
            SimulationModel.CompliantJointActuator[],
    )
end

"""Return whether any configured effector carries a robot-arm plan at all."""
function _any_robot_arm_effector(args)::Bool
    if hasproperty(args, :control_model) && hasproperty(args.control_model, :control_effectors)
        @inbounds for effector in args.control_model.control_effectors
            hasproperty(effector, :plan) &&
                getproperty(effector, :plan) isa SimulationModel.RobotArmPlan && return true
        end
    end
    if hasproperty(args, :dynamics_model) && hasproperty(args.dynamics_model, :dynamic_effectors)
        @inbounds for effector in args.dynamics_model.dynamic_effectors
            hasproperty(effector, :plan) &&
                getproperty(effector, :plan) isa SimulationModel.RobotArmPlan && return true
        end
    end
    return false
end

"""Cached robot-arm presence check for the RHS hot path (lazy, run-constant)."""
@inline function _robot_arm_present(p)::Bool
    (p !== nothing && hasproperty(p, :shared_buffers) &&
        hasproperty(p.shared_buffers, :robot_arm_present)) || return true
    cached = p.shared_buffers.robot_arm_present[]
    cached === nothing || return cached
    present = _any_robot_arm_effector(p.args)
    p.shared_buffers.robot_arm_present[] = present
    return present
end

"""Find the active robot-arm coupling configuration for one spacecraft."""
function _robot_arm_coupling(args, sat_idx::Int, t::Float64)
    if hasproperty(args, :control_model) && hasproperty(args.control_model, :control_effectors)
        @inbounds for effector in args.control_model.control_effectors
            _robot_arm_effector_matches(effector, sat_idx) && return _robot_arm_coupling_from_effector(effector, t)
        end
    end
    if hasproperty(args, :dynamics_model) && hasproperty(args.dynamics_model, :dynamic_effectors)
        @inbounds for effector in args.dynamics_model.dynamic_effectors
            _robot_arm_effector_matches(effector, sat_idx) && return _robot_arm_coupling_from_effector(effector, t)
        end
    end
    return nothing
end

"""Apply coupled cloth robot-arm state derivatives to one spacecraft RHS view."""
@inline function _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, sat_idx::Int, t::Float64, forces, torques)
    # Fast path: skip the per-satellite effector-tuple scan when the run has no
    # robot-arm effector anywhere (cached in shared_buffers on first use).
    _robot_arm_present(p) || return nothing
    coupling = _robot_arm_coupling(p.args, sat_idx, t)
    coupling === nothing && return nothing
    SimulationModel.assign_coupled_cloth_robot_arm_rhs!(
        du_view,
        sc_view,
        coupling.plan,
        coupling.t_s,
        forces,
        torques;
        k_translation_n_m=coupling.k_translation_n_m,
        c_translation_n_s_m=coupling.c_translation_n_s_m,
        k_rotation_n_m_rad=coupling.k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=coupling.c_rotation_n_m_s_rad,
        joint_actuators=coupling.joint_actuators,
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
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(p, length(spacecraft))
    if use_rhs_batch
        _prefill_shared_body_samples!(p, t, sc_state, dynamic_effectors)
    end
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
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(p, length(spacecraft))
    if use_rhs_batch
        _prefill_shared_body_samples!(p, t, sc_state, dynamic_effectors)
    end
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
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
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
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics_slow!

@inline function _drag_state_buffer_current(p, sat_idx::Int, t::Float64)::Bool
    times = p.shared_buffers.in_atmosphere_sample_t
    return sat_idx <= length(times) && times[sat_idx] == t
end

@inline function _spacecraft_outside_atmosphere_for_current_state(sc, p, sat_idx::Int, t::Float64)::Bool
    # The above-EI shortcut is only sound when the density model is guaranteed
    # to be vacuum there. For every other model EI is a step-size/tolerance
    # boundary, not a force gate: an orbit that never crosses EI downward must
    # still see continuous aero (negligible densities zero out inside the
    # wrench itself via its rho <= eps short-circuit).
    density_model = SimulationModel.SimulationCallbacks._density_model_for_sat(p, sat_idx)
    SimulationModel.EnvironmentModels.density_vanishes_above_entry_interface(density_model) || return false

    if _drag_state_buffer_current(p, sat_idx, t)
        return !p.shared_buffers.in_atmosphere[sat_idx]
    end

    planet_frame = p.shared_buffers.rhs_planet_frame_prefilled[] ?
        sample_buffered_planet_frame(p, sat_idx) :
        sample_planet_frame(sc, p, sat_idx, t)
    isfinite(planet_frame.alt_m) || return false
    return planet_frame.alt_m > p.args.environment_model.EI * 1e3
end

@inline function _all_active_spacecraft_outside_atmosphere(sc_state, p, t::Float64)::Bool
    @inbounds for i in eachindex(sc_state)
        p.is_active[i] || continue
        _spacecraft_outside_atmosphere_for_current_state(sc_state[i], p, i, t) || return false
    end
    return true
end

function spacecraft_dynamics_implicit_atmosphere!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    p.shared_buffers.current_time[] = t
    # Fast path: if the current state proves that no active spacecraft is inside
    # the atmosphere, all implicit drag terms are zero.  We only trust the staged
    # drag-state flag when it was stamped for this RHS time; otherwise we derive
    # the decision from the current altitude to avoid skipping first-step or
    # direct RHS evaluations that start inside the entry interface.
    if _all_active_spacecraft_outside_atmosphere(sc_state, p, t)
        du .= 0.0
        return nothing
    end
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
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(p, length(spacecraft))
    if use_rhs_batch
        _prefill_shared_body_samples!(p, t, sc_state, dynamic_effectors)
    end
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i] || _spacecraft_outside_atmosphere_for_current_state(sc_state[i], p, i, t)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i] || _spacecraft_outside_atmosphere_for_current_state(sc_state[i], p, i, t)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
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
    use_rhs_batch = plan.mode != :serial && _rhs_batch_parallel_enabled(p, length(spacecraft))
    if use_rhs_batch
        _prefill_shared_body_samples!(p, t, sc_state, dynamic_effectors)
    end
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
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
    use_rhs_batch = _rhs_batch_parallel_enabled(p, length(spacecraft))
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
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)

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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
                rw_torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, rw_torque_body, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)

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
                        rw_assembly=spacecraft[i].root.rw_assembly,
                        rw_torque_body=SVector{3, Float64}(rw_torque_body),
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
    sc_shapes = map(eachindex(args.dynamics_model.spacecraft)) do i
        sc = args.dynamics_model.spacecraft[i]
        n_bodies = length(sc.links)
        mass = sc.dry_mass + sc.prop_mass
        base_shape = if args.mission_configuration.orientation_sim
            (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies),
                q = Float64[0.0, 0.0, 0.0, 1.0], 
                ω = zeros(3)
            )
        else
            (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies)
            )
        end
        n_rw = args.mission_configuration.orientation_sim ? sc.root.rw_assembly.n_wheels : 0
        base_shape = n_rw > 0 ? merge(base_shape, (h_wheels = zeros(n_rw),)) : base_shape
        coupling = _robot_arm_coupling(args, i, 0.0)
        coupling === nothing && return base_shape
        return merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.plan))
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
            if spacecraft.root.rw_assembly.n_wheels > 0
                sc_view.h_wheels .= spacecraft.root.rw_assembly.h_wheels
            end
        end
        coupling = _robot_arm_coupling(args, i, 0.0)
        if coupling !== nothing
            SimulationModel.initialize_coupled_cloth_robot_arm_state!(sc_view, coupling.plan; t_s=coupling.t_s)
        end
    end

    return state
end
