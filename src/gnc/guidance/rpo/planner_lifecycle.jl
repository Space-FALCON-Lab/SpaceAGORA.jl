"""
    RPOPlanningEvent(time_s; reason=:forced, token=1)

Deliver a request at the first accepted guidance tick at or after `time_s`.
Reasons are `:forced`, `:replan`, or `:retime`. Positive tokens identify events;
repeated delivery of one token is idempotent. Distinct tokens at one tick are
separate requests. Retiming must be supported explicitly by the planner.
"""
struct RPOPlanningEvent
    time_s::Float64
    reason::Symbol
    token::Int
    function RPOPlanningEvent(time_s; reason::Symbol=:forced, token::Integer=1)
        isfinite(time_s) && time_s >= 0 || throw(ArgumentError("Event time must be finite and nonnegative."))
        reason in (:forced, :replan, :retime) || throw(ArgumentError("Unsupported planning event."))
        token > 0 || throw(ArgumentError("Event token must be positive."))
        new(Float64(time_s), reason, Int(token))
    end
end

"""Planning failure with copied run diagnostics, original cause and backtrace."""
struct RPOPlanningError <: Exception
    chaser_id::Int
    request_id::Int
    reason::Symbol
    report::Any
    cause::Any
    backtrace::Any
end
function Base.showerror(io::IO, e::RPOPlanningError)
    print(io, "RPO planning failed for spacecraft ", e.chaser_id,
        ", request ", e.request_id, ": ", e.reason)
    e.cause === nothing || (print(io, ". Cause: "); showerror(io, e.cause))
end

mutable struct PlannerRuntime
    state::Any
    rng::MersenneTwister
    stream_words::Vector{UInt32}
    request_id::Int
    installed_request::Any
    reference::Union{Nothing,P.RPOReference}
    processed_tokens::Set{Int}
    last_guidance_time::Float64
    last_planning_time::Float64
    records::Vector{NamedTuple}
end

# Kept separate from RPOGuidanceModel and RPOPlanBuffer for legacy compatibility.
Base.@kwdef mutable struct PlannerGuidance <: S.AbstractGuidanceModel
    planner::P.AbstractRPOPlanner
    chaser_id::Int
    target_id::Int
    chaser_idx::Int
    target_idx::Int
    goal_rtn_m::NTuple{3,Float64}
    geometry::S.RPOReferenceGeometry
    geometry_revision::String
    constraints::P.RPOPlanningConstraints
    validation::P.RPOValidationSettings
    reference_dt_s::Float64
    preview_horizon_steps::Int
    validity_s::Float64
    seed::Int
    events::Vector{RPOPlanningEvent}
    replan_interval_s::Float64
    tracking_error_limit_m::Float64
    plan_buffer::S.RPOPlanBuffer
    runtime::Union{Nothing,PlannerRuntime} = nothing
end

# SHA-256 over this versioned ASCII encoding, decoded little-endian, not hash().
function stream_words(seed, chaser_id, target_id)
    digest = sha256("SpaceAGORA/RPO/v1/$(seed)/$(chaser_id)/$(target_id)")
    return [sum(UInt32(digest[i+j]) << (8j) for j in 0:3) for i in 1:4:32]
end

function report_binding(g::PlannerGuidance)
    r = g.runtime
    return deepcopy((chaser_id=g.chaser_id, target_id=g.target_id,
        planner=string(typeof(g.planner)), settings=g.planner, seed=g.seed,
        generator="MersenneTwister", stream_mapping="SHA256/SpaceAGORA/RPO/v1/little-endian",
        stream_words=r === nothing ? stream_words(g.seed,g.chaser_id,g.target_id) : r.stream_words,
        state_source=:truth, frame=:target_rtn, geometry_revision=g.geometry_revision,
        constraints=g.constraints, validation_settings=g.validation,
        request_id=r === nothing ? 0 : r.request_id,
        active_reference=r === nothing ? nothing : r.reference,
        records=r === nothing ? NamedTuple[] : r.records))
end
function fail!(g, reason; cause=nothing, backtrace=nothing, details=NamedTuple())
    r = g.runtime
    r === nothing || push!(r.records, (event=:failure, request_id=r.request_id,
        reason=reason, details=deepcopy(details)))
    throw(RPOPlanningError(g.chaser_id, r === nothing ? 0 : r.request_id,
        reason, report_binding(g), cause, backtrace))
end

function L.preflight_guidance(g::PlannerGuidance, args; isolate_state)
    isolate_state || throw(ArgumentError("The RPO planner lifecycle requires isolate_state=true."))
    ss = args.simulation_settings
    (ss.checkpoint_enabled || ss.resume_from_checkpoint) &&
        throw(ArgumentError("The RPO planner lifecycle does not support checkpoint writing or resume."))
    ids = [sc.id for sc in args.dynamics_model.spacecraft]
    length(ids) == length(unique(ids)) || throw(ArgumentError("Spacecraft IDs must be unique."))
    g.chaser_id in ids && g.target_id in ids || throw(ArgumentError("Planner spacecraft ID is missing."))
    count(x -> x isa PlannerGuidance && x.chaser_id == g.chaser_id,
        args.guidance_model.guidance_effectors) == 1 || throw(ArgumentError("One planner owner per chaser is required."))
    controls = [c for c in args.control_model.control_effectors if
        c isa S.RPOMPCControlModel && c.plan_buffer === g.plan_buffer]
    length(controls) == 1 || throw(ArgumentError("Planner and one controller must share their reference buffer."))
    c = only(controls)
    c.control_dt_s == g.reference_dt_s && c.controller.horizon == g.preview_horizon_steps ||
        throw(ArgumentError("Planner and controller preview settings disagree."))
    return nothing
end

function L.initialize_guidance!(g::PlannerGuidance, u, p, t)
    ids = [sc.id for sc in p.args.dynamics_model.spacecraft]
    g.chaser_idx = only(findall(==(g.chaser_id), ids))
    g.target_idx = only(findall(==(g.target_id), ids))
    c = only(c for c in p.args.control_model.control_effectors if
        c isa S.RPOMPCControlModel && c.plan_buffer === g.plan_buffer)
    c.chaser_idx, c.target_idx = g.chaser_idx, g.target_idx
    c.held = S.RPOHeldActuation()
    c.command_log = S.RPOControlCommandLog()
    fill!(c.controller.U_prev, 0)
    words = stream_words(g.seed,g.chaser_id,g.target_id)
    g.runtime = PlannerRuntime(nothing, MersenneTwister(words), words, 0,
        nothing, nothing, Set{Int}(), -Inf, -Inf, NamedTuple[])
    g.plan_buffer.valid = false
    try
        g.runtime.state = P.initialize_planner(g.planner,
            (seed=g.seed, chaser_id=g.chaser_id, target_id=g.target_id,
             epoch=deepcopy(p.args.initial_time), stream_words=copy(words)))
        request_plan!(g, u, p, Float64(t), :initial)
    catch e
        e isa RPOPlanningError && rethrow()
        fail!(g, :initialization_exception; cause=e, backtrace=catch_backtrace())
    end
    return nothing
end

function request_plan!(g, u, p, t, reason)
    r = g.runtime::PlannerRuntime
    r.request_id += 1
    try
        a, b = u.sc[g.chaser_idx], u.sc[g.target_idx]
        x = S.FrameTransforms.inertial_to_rtn_relative_state(a.pos,a.vel,b.pos,b.vel)
        # Retain a trusted snapshot. Planner code only receives another copy.
        req = P.RPOPlanningRequest(request_id=r.request_id, reason=reason,
            chaser_id=g.chaser_id, target_id=g.target_id, epoch=p.args.initial_time,
            time_s=t, x_rtn=x, target_state_ii=vcat(b.pos,b.vel),
            goal_rtn_m=g.goal_rtn_m, geometry=g.geometry,
            geometry_revision=g.geometry_revision, constraints=g.constraints,
            reference_dt_s=g.reference_dt_s, preview_horizon_steps=g.preview_horizon_steps,
            valid_until_s=t+g.validity_s, validation=g.validation)
        push!(r.records, (event=:request, time_s=t, request=deepcopy(req)))
        caps = P.validate_rpo_capabilities(g.planner, req)
        caps.accepted || fail!(g, caps.reason)
        result = reason === :retime ?
            P.retime_rpo!(r.state, g.planner, deepcopy(req), deepcopy(r.reference), r.rng) :
            P.plan_rpo!(r.state, g.planner, deepcopy(req), r.rng)
        result isa P.RPOPlanningResult || fail!(g, :invalid_result_type)
        result = deepcopy(result)
        checked = P.validate_rpo_result(req, result;
            clearance_at=(point,geom)->S.rpo_clearance_distance_to_station(point,geom))
        checked.accepted || fail!(g, checked.reason;
            details=(termination=result.termination, diagnostics=result.diagnostics, validation=checked))
        ref = deepcopy(result.reference)
        # Diagnostics never choose the controller's geometry. This is the validated polyline.
        plan = S.RPOPlan(valid=true, t_ref_s=copy(ref.t_ref_s),
            r_ref_rtn=copy(ref.r_ref_rtn_m), v_ref_rtn=copy(ref.v_ref_rtn_mps),
            path_rtn=copy(ref.r_ref_rtn_m), cost=NaN,
            diagnostics=(request_id=req.request_id, validation=checked,
                         planner=deepcopy(result.diagnostics)))
        record = deepcopy((event=:installed, time_s=t, request_id=req.request_id,
            reason=reason, termination=result.termination, validation=checked,
            reference=ref, diagnostics=result.diagnostics))
        # No yield or user callback between preparing and committing the complete plan.
        S.update_rpo_plan_buffer!(g.plan_buffer, plan, t)
        r.reference, r.installed_request, r.last_planning_time = ref, req, t
        push!(r.records, record)
    catch e
        e isa RPOPlanningError && rethrow()
        fail!(g, :planner_exception; cause=e, backtrace=catch_backtrace())
    end
    return nothing
end

function S.GuidanceHooks.calcGuidanceEffect!(g::PlannerGuidance, u, p, t::Float64, sat_idx::Int)
    sat_idx == g.chaser_idx || return nothing
    r = g.runtime
    r === nothing && fail!(g, :uninitialized)
    # Same accepted tick may arrive through periodic and thruster callbacks.
    duplicate_tick = t == r.last_guidance_time
    due = any(e -> e.time_s <= t + g.validation.time_atol_s &&
        !(e.token in r.processed_tokens), g.events)
    duplicate_tick && !due && return nothing
    t >= r.last_guidance_time || fail!(g, :nonmonotonic_guidance_time)
    r.last_guidance_time = t
    push!(r.records, (event=:guidance, time_s=t, request_id=r.request_id))
    requested = false
    for e in g.events
        e.time_s <= t + g.validation.time_atol_s && !(e.token in r.processed_tokens) || continue
        push!(r.processed_tokens, e.token)
        request_plan!(g,u,p,t,e.reason)
        requested = true
    end
    if !requested && !duplicate_tick
        a,b = u.sc[g.chaser_idx],u.sc[g.target_idx]
        x = S.FrameTransforms.inertial_to_rtn_relative_state(a.pos,a.vel,b.pos,b.vel)
        preview = S.ControlHooks.rpo_ref_preview(g.plan_buffer.plan,
            max(0.0,t-g.plan_buffer.updated_at_s),g.reference_dt_s,1)
        if t-r.last_planning_time >= g.replan_interval_s ||
                norm(x[1:3]-preview[1:3,1]) > g.tracking_error_limit_m
            request_plan!(g,u,p,t,:replan)
        end
    end
    return nothing
end

function L.before_reference_control!(g::PlannerGuidance, c, u, p, t)
    c isa S.RPOMPCControlModel && c.plan_buffer === g.plan_buffer || return nothing
    r = g.runtime
    r === nothing && fail!(g,:uninitialized)
    r.reference === nothing && fail!(g,:missing_reference)
    P.rpo_reference_is_current(r.installed_request,r.reference,t) ||
        fail!(g,:expired_reference; details=(time_s=t,valid_until_s=r.reference.valid_until_s,
            preview_end_s=t+g.preview_horizon_steps*g.reference_dt_s))
    push!(r.records,(event=:control, time_s=Float64(t), request_id=r.request_id))
    return nothing
end

"""
    rpo_run_report(solution)

Return copied planner requests, validation, references, stream identities,
controller commands and sampled states from an isolated `run_simulation` result.
This pilot reports reference acceptance; it does not certify physical tracking.
"""
function rpo_run_report(solution)
    args = solution.prob.p.args
    bindings = [g for g in args.guidance_model.guidance_effectors if g isa PlannerGuidance]
    return deepcopy((planners=[report_binding(g) for g in bindings],
        commands=[(chaser_id=g.chaser_id, log=only(c.command_log for c in
            args.control_model.control_effectors if c isa S.RPOMPCControlModel &&
            c.plan_buffer === g.plan_buffer)) for g in bindings],
        times_s=solution.t, states=solution.u))
end

# Analytic clearance of a segment to the same point-sphere approximation used
# by the trusted sample validator. This is not a watertight CAD collision model.
function station_segment_clearance(a,b,geometry)
    d = collect(b)-collect(a); denom = dot(d,d)
    return minimum(begin
        q = denom == 0 ? 0.0 : clamp(dot(point-collect(a),d)/denom,0.0,1.0)
        norm(collect(a)+q*d-point)
    end for point in eachcol(geometry.station.points_body)) -
        geometry.station.keepout_radius_m - maximum(geometry.chaser.half_extents_body)
end
