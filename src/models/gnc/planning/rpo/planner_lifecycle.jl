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

"""Planning failure with copied run diagnostics, original cause and backtrace; a `GNCAlgorithms.GNCAlgorithmError`."""
struct RPOPlanningError <: GA.GNCAlgorithmError
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

# The planner's run-owned state is the shared lifecycle record; `product` holds
# the installed reference.
const PlannerRuntime = GA.AlgorithmRuntime

# Run-owned trigger state: whether the installed reference is the family's hold
# (`:none`, `:holding`, or `:retried` once a retry has been rejected) and the
# geometry revision the installed reference was last checked against.
mutable struct TriggerState
    hold::Symbol
    revision::String
end

# The fields after `runtime` configure the contract extensions; their defaults are
# the static geometry, the `:fail` policy and no retiming trigger.
Base.@kwdef mutable struct PlannerGuidance <: S.AbstractGuidanceModel
    planner::P.AbstractRPOPlanner
    chaser_id::Int
    target_id::Int
    chaser_index::Int
    target_index::Int
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
    family::P.RPOReferenceFamily = P.RPOReferenceFamily()
    geometry_provider::Union{Nothing,P.AbstractRPOGeometryProvider} = nothing
    retime_tracking_error_m::Float64 = Inf
    hold_position_tolerance_m::Float64 = .02
    hold_speed_tolerance_mps::Float64 = .002
    triggers::Union{Nothing,TriggerState} = nothing
end

# The planner ticks for its chaser and reads the target (GNCRuntime). The indices
# are resolved from the spacecraft IDs at the run's start.
GR.spacecraft_scope(g::PlannerGuidance) = GR.CoupledSpacecraft(g.chaser_index, (g.target_index,))

# SHA-256 over "SpaceAGORA/RPO/v1/<seed>/<chaser_id>/<target_id>", decoded
# little-endian, not hash(): the shared stream mapping with the RPO namespace.
stream_words(seed, chaser_id, target_id) = GA.gnc_stream_words("RPO/v1", seed, chaser_id, target_id)

function report_binding(g::PlannerGuidance)
    r = g.runtime
    return deepcopy((chaser_id=g.chaser_id, target_id=g.target_id,
        planner=string(typeof(g.planner)), settings=g.planner, seed=g.seed,
        generator="MersenneTwister", stream_mapping="SHA256/SpaceAGORA/RPO/v1/little-endian",
        stream_words=r === nothing ? stream_words(g.seed,g.chaser_id,g.target_id) : r.stream_words,
        state_source=:truth, frame=:target_rtn, geometry_revision=g.geometry_revision,
        constraints=g.constraints, validation_settings=g.validation,
        request_id=r === nothing ? 0 : r.request_id,
        active_reference=r === nothing ? nothing : r.product,
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
    P.require_planner_support(g.planner)
    isolate_state || throw(ArgumentError("The RPO planner lifecycle requires isolate_state=true."))
    ids = [sc.id for sc in args.dynamics_model.spacecraft]
    length(ids) == length(unique(ids)) || throw(ArgumentError("Spacecraft IDs must be unique."))
    g.chaser_id in ids && g.target_id in ids || throw(ArgumentError("Planner spacecraft ID is missing."))
    count(x -> x isa PlannerGuidance && x.chaser_id == g.chaser_id,
        args.guidance_model.guidance_effectors) == 1 || throw(ArgumentError("One planner owner per chaser is required."))
    # Its pairing with the controller is a check of the plan-buffer link (below).
    return nothing
end

# The planner writes the plan buffer that its controller reads (GNCRuntime link).
GR.gnc_links(g::PlannerGuidance) = (GR.Produces(g.plan_buffer),)

# A plan buffer written by a planner: exactly one RPO MPC controller reads it, with
# the planner's preview settings. A buffer no planner writes gets no check.
function GR.check_link(b::S.RPOPlanBuffer, producers, consumers, args)
    planners = [x for x in producers if x isa PlannerGuidance]
    isempty(planners) && return nothing
    controls = [c for c in consumers if c isa S.RPOMPCControlModel]
    length(controls) == 1 || throw(ArgumentError("Planner and one controller must share their reference buffer."))
    c = only(controls)
    for g in planners
        c.control_dt_s == g.reference_dt_s && c.controller.horizon == g.preview_horizon_steps ||
            throw(ArgumentError("Planner and controller preview settings disagree."))
    end
    return nothing
end

# At the run's start the planner's controller is bound to the planner's spacecraft
# and its run state reset, before either model's own initializer.
function GR.link_initialize!(b::S.RPOPlanBuffer, producers, consumers, u, p, t)
    planners = [x for x in producers if x isa PlannerGuidance]
    isempty(planners) && return nothing
    g = only(planners)
    c = only(c for c in consumers if c isa S.RPOMPCControlModel)
    ids = [sc.id for sc in p.args.dynamics_model.spacecraft]
    c.chaser_index = only(findall(==(g.chaser_id), ids))
    c.target_index = only(findall(==(g.target_id), ids))
    c.held = S.RPOHeldActuation()
    c.command_log = S.RPOControlCommandLog()
    fill!(c.controller.U_prev, 0)
    return nothing
end

function L.initialize_guidance!(g::PlannerGuidance, u, p, t)
    ids = [sc.id for sc in p.args.dynamics_model.spacecraft]
    g.chaser_index = only(findall(==(g.chaser_id), ids))
    g.target_index = only(findall(==(g.target_id), ids))
    words = stream_words(g.seed,g.chaser_id,g.target_id)
    g.runtime = GA.AlgorithmRuntime(words)
    g.triggers = TriggerState(:none, g.geometry_revision)
    g.plan_buffer.valid = false
    # A resumed run puts back the planner's run state from its checkpoint right
    # after the run-start initializers (`restore_checkpoint_state!`): no initial request.
    hasproperty(p, :shared_buffers) && p.shared_buffers.model_resume[] && return nothing
    _guarded_initialization(g) do
        g.runtime.state = GA.initialize_algorithm(g.planner,
            (seed=g.seed, chaser_id=g.chaser_id, target_id=g.target_id,
             epoch=deepcopy(p.args.initial_time), stream_words=copy(words)))
    end
    # A checkpoint holds the algorithm's own state only for a planner that declares
    # restart. A run set to resume writes checkpoints too.
    ss = p.args.simulation_settings
    g.runtime.state === nothing || !(ss.checkpoint_enabled || ss.resume_from_checkpoint) ||
        P.planner_capabilities(g.planner).restart || throw(ArgumentError(
            "$(nameof(typeof(g.planner))) keeps run state of its own (initialize_planner) and does not declare " *
            "restart in planner_capabilities, so a checkpoint cannot hold it; run without checkpoints."))
    _guarded_initialization(() -> request_plan!(g, u, p, Float64(t), :initial), g)
    return nothing
end

# A failure while the planner initializes is the run's planning failure.
function _guarded_initialization(f, g)
    try
        f()
    catch e
        e isa RPOPlanningError && rethrow()
        fail!(g, :initialization_exception; cause=e, backtrace=catch_backtrace())
    end
    return nothing
end

# The planner's run state goes into checkpoint files: the request driver's runtime
# (the random stream, the installed request and reference, the records and the
# algorithm's own state), the trigger state and the plan buffer the planner shares
# with its controller.
function L.checkpoint_state(g::PlannerGuidance)
    g.runtime === nothing && return nothing
    b = g.plan_buffer
    return deepcopy((runtime=g.runtime, triggers=g.triggers,
        plan_buffer=(valid=b.valid, plan=b.plan, updated_at_s=b.updated_at_s)))
end
function L.restore_checkpoint_state!(g::PlannerGuidance, data)
    g.runtime, g.triggers = data.runtime, data.triggers
    b, saved = g.plan_buffer, data.plan_buffer
    b.valid, b.plan, b.updated_at_s = saved.valid, saved.plan, saved.updated_at_s
    return nothing
end

# The planner is the RPO binding of the shared request driver (GNCAlgorithms):
# the driver runs request, capability check, proposal, validation and install in
# a fixed order; the methods below supply the RPO pieces.
GA.algorithm_runtime(g::PlannerGuidance) = g.runtime
GA.bound_algorithm(g::PlannerGuidance) = g.planner
GA.bound_family(g::PlannerGuidance) = g.family
GA.bound_link(g::PlannerGuidance) = g.plan_buffer
GA.trusted_queries(::PlannerGuidance) =
    (clearance_at=(point,geom)->S.rpo_clearance_distance_to_station(point,geom),)
GA.fail_algorithm_request!(g::PlannerGuidance, reason; kwargs...) = fail!(g, reason; kwargs...)
function GA.build_algorithm_request(g::PlannerGuidance, u, p, t, reason, request_id)
    a, b = u.sc[g.chaser_index], u.sc[g.target_index]
    x = S.FrameTransforms.inertial_to_rtn_relative_state(a.pos,a.vel,b.pos,b.vel)
    # A provider supplies this request's owned geometry snapshot and its revision.
    provider = g.geometry_provider
    geometry, revision = provider === nothing ? (g.geometry, g.geometry_revision) :
        (P.rpo_geometry_snapshot(provider, g.geometry, t), P.rpo_geometry_revision(provider, g.geometry_revision, t))
    return P.RPOPlanningRequest(request_id=request_id, reason=reason,
        chaser_id=g.chaser_id, target_id=g.target_id, epoch=p.args.initial_time,
        time_s=t, x_rtn=x, target_state_ii=vcat(b.pos,b.vel),
        goal_rtn_m=g.goal_rtn_m, geometry=geometry,
        geometry_revision=revision, constraints=g.constraints,
        reference_dt_s=g.reference_dt_s, preview_horizon_steps=g.preview_horizon_steps,
        valid_until_s=t+g.validity_s, validation=g.validation)
end

# Install the validated reference: the buffer the planner shares with its controller.
function GA.install!(buffer::S.RPOPlanBuffer, ::P.RPOReferenceFamily, ref::P.RPOReference,
        request::P.RPOPlanningRequest, t; validation, diagnostics)
    # Diagnostics never choose the controller's geometry. This is the validated polyline.
    plan = S.RPOPlan(valid=true, t_ref_s=copy(ref.t_ref_s),
        r_ref_rtn=copy(ref.r_ref_rtn_m), v_ref_rtn=copy(ref.v_ref_rtn_mps),
        path_rtn=copy(ref.r_ref_rtn_m), cost=NaN,
        diagnostics=(request_id=request.request_id, validation=validation,
                     planner=deepcopy(diagnostics)))
    S.update_rpo_plan_buffer!(buffer, plan, t)
    return nothing
end

# `hold_retry` marks a hold's retry or renewal: if the family installs another hold,
# it waits for its renewal, a new geometry revision or an event, not for settling.
function request_plan!(g, u, p, t, reason; hold_retry::Bool=false)
    GA.run_algorithm_request!(g, u, p, t, reason; exception_reason=:planner_exception)
    r, s = g.runtime, g.triggers
    held = last(r.records).event === :fallback_installed
    s.hold = !held ? :none : hold_retry ? :retried : :holding
    s.revision = r.installed_request.geometry_revision
    return nothing
end

trigger_record(t, trigger; details...) = (event=:trigger, time_s=t, trigger=trigger, details=(; details...))

function provider_call(g, f, base, t)
    try
        return f(g.geometry_provider, base, t)
    catch e
        fail!(g, :geometry_provider_exception; cause=e, backtrace=catch_backtrace())
    end
end

# The smallest trusted clearance along the installed reference's remaining polyline,
# from its position at `t` to its end, sampled as validation samples a reference.
# `-Inf` when the samples would exceed the validation work cap.
function remaining_clearance(g, geometry, t)
    ref, settings = g.runtime.product, g.validation
    times, positions = ref.t_ref_s, ref.r_ref_rtn_m
    n = length(times); elapsed = t - ref.origin_time_s
    j = clamp(searchsortedlast(times, elapsed), 1, n)
    start = j == n ? positions[:,n] :
        (q = clamp((elapsed-times[j])/(times[j+1]-times[j]), 0., 1.); (1-q) .* positions[:,j] .+ q .* positions[:,j+1])
    points = hcat(start, positions[:,j+1:n])
    clearance_at = GA.trusted_queries(g).clearance_at
    smallest = Float64(clearance_at(Tuple(points[:,1]), geometry)); count = 1
    for i in 2:size(points,2)
        a, b = points[:,i-1], points[:,i]
        steps = norm(b-a)/settings.clearance_sample_ds_m
        isfinite(steps) && steps <= settings.max_clearance_samples - count || return -Inf
        steps = max(1, ceil(Int, steps)); count += steps
        count <= settings.max_clearance_samples || return -Inf
        for k in 1:steps
            smallest = min(smallest, Float64(clearance_at(Tuple((1-k/steps) .* a .+ (k/steps) .* b), geometry)))
        end
    end
    return smallest
end

# A new geometry revision replans a hold, and replans another reference only if its
# remaining path no longer clears the new geometry. Returns whether it requested.
function revision_requested!(g, u, p, t)
    g.geometry_provider === nothing && return false
    r, s = g.runtime, g.triggers
    revision = provider_call(g, P.rpo_geometry_revision, g.geometry_revision, t)
    revision isa AbstractString && !isempty(revision) || fail!(g, :invalid_geometry_revision)
    revision == s.revision && return false
    smallest = NaN
    if s.hold === :none
        smallest = remaining_clearance(g, provider_call(g, P.rpo_geometry_snapshot, g.geometry, t), t)
        if smallest >= g.constraints.clearance_m
            s.revision = revision
            push!(r.records, trigger_record(t, :geometry_revision; revision=String(revision),
                action=:kept, min_clearance_m=smallest))
            return false
        end
    end
    push!(r.records, trigger_record(t, :geometry_revision; revision=String(revision),
        action=:replan, min_clearance_m=smallest))
    request_plan!(g,u,p,t,:replan)
    return true
end

# Requests beyond the interval and tracking-error replans, at most one per tick:
# - a new geometry revision (`revision_requested!`);
# - a hold: renewed before it expires, and retried once the chaser settles on it;
# - otherwise, a tracking error above `retime_tracking_error_m` asks for a retime.
function request_triggered!(g, u, p, t, x, preview)
    revision_requested!(g, u, p, t) && return nothing
    r, s = g.runtime, g.triggers
    if s.hold !== :none
        if !GA.is_current(g.family, r.installed_request, r.product, t + g.reference_dt_s)
            push!(r.records, trigger_record(t, :hold_renewal))
            return request_plan!(g,u,p,t,:replan; hold_retry=true)
        elseif s.hold === :holding && norm(x[1:3]-r.product.r_ref_rtn_m[:,1]) <= g.hold_position_tolerance_m &&
                norm(x[4:6]) <= g.hold_speed_tolerance_mps
            push!(r.records, trigger_record(t, :hold_settled))
            return request_plan!(g,u,p,t,:replan; hold_retry=true)
        end
    elseif norm(x[1:3]-preview[1:3,1]) > g.retime_tracking_error_m
        push!(r.records, trigger_record(t, :retime_tracking_error))
        return request_plan!(g,u,p,t,:retime)
    end
    return nothing
end

function S.GuidanceHooks.update_guidance!(g::PlannerGuidance, u, p, t::Float64, sat_idx::Int)
    sat_idx == g.chaser_index || return nothing
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
        a,b = u.sc[g.chaser_index],u.sc[g.target_index]
        x = S.FrameTransforms.inertial_to_rtn_relative_state(a.pos,a.vel,b.pos,b.vel)
        preview = S.ControlHooks.rpo_ref_preview(g.plan_buffer.plan,
            max(0.0,t-g.plan_buffer.updated_at_s),g.reference_dt_s,1)
        # Guidance tick subtraction can round just below the requested interval.
        # Use the same time tolerance as scheduled events, without widening expiry.
        # While a hold is installed, interval and tracking-error replans wait: braking
        # onto it is expected, and only the hold's own triggers replace it, so a
        # planner that keeps refusing is asked again at a bounded rate.
        if g.triggers.hold === :none && (t-r.last_request_time >= g.replan_interval_s - g.validation.time_atol_s ||
                norm(x[1:3]-preview[1:3,1]) > g.tracking_error_limit_m)
            request_plan!(g,u,p,t,:replan)
        else
            request_triggered!(g,u,p,t,x,preview)
        end
    end
    return nothing
end

function L.before_reference_control!(g::PlannerGuidance, c, u, p, t)
    c isa S.RPOMPCControlModel && c.plan_buffer === g.plan_buffer || return nothing
    r = g.runtime
    r === nothing && fail!(g,:uninitialized)
    r.product === nothing && fail!(g,:missing_reference)
    GA.is_current(g.family,r.installed_request,r.product,t) ||
        fail!(g,:expired_reference; details=(time_s=t,valid_until_s=r.product.valid_until_s,
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
