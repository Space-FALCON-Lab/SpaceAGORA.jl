# Experimental contact-load orchestration. No provider/native code belongs here.
const CLS = ContactLoadServices

function CLS.with_contact_model(args::SimulationConfiguration, model::CLS.AbstractContactModel)
    return SimulationModel.SimConfig._with_configuration(args; contact_models=(model,))
end

function _contact_bodies(args)
    sc = only(args.dynamics_model.spacecraft)
    body = only(only(sc.attachments).model.bodies)
    return (CLS.ContactBody(sc.dry_mass + sc.prop_mass, sc.inertia_tensor),
        CLS.ContactBody(body.mass_kg, body.inertia_body_kg_m2))
end

function _validate_contact_coupling!(args, solver_mode)
    isempty(args.contact_models) && return nothing
    length(args.contact_models) == 1 && only(args.contact_models) isa CLS.AbstractContactModel ||
        throw(ArgumentError("Contact coupling requires one AbstractContactModel."))
    length(args.dynamics_model.spacecraft) == 1 || throw(ArgumentError("Contact coupling supports one bus and one attachment body."))
    solver_mode === :tsit5 && !_run_solver_config(args).parallel ||
        throw(ArgumentError("Contact coupling requires serial :tsit5 integration."))
    isempty(args.external_propagators) && args.interlink_model === nothing ||
        throw(ArgumentError("Contact coupling cannot be combined with external propagation or interlinks."))
    isempty(args.dynamics_model.dynamic_effectors) && isempty(args.control_model.control_effectors) &&
        isempty(args.guidance_model.guidance_effectors) && isempty(args.navigation_model.navigation_effectors) ||
        throw(ArgumentError("Initial contact coupling supports gravity-free, unactuated runs without other effectors."))
    settings = args.simulation_settings
    !settings.checkpoint_enabled && !settings.resume_from_checkpoint ||
        throw(ArgumentError("Contact coupling does not support checkpoint/resume."))
    args.mission_configuration.mission_type == SimulationModel.MissionTime ||
        throw(ArgumentError("Contact coupling requires a time-bounded mission."))
    tol = args.integration_tolerances
    all(h -> isfinite(h) && 0 < h <= 2.0^-10, (tol.dt_max_orbit, tol.dt_max_atmosphere)) ||
        throw(ArgumentError("Contact coupling requires both maximum step caps <= 2^-10 s; start with 2^-12 s."))
    sc = only(args.dynamics_model.spacecraft)
    !_spacecraft_articulated(sc) && length(sc.attachments) == 1 && length(sc.links) == 1 ||
        throw(ArgumentError("Contact coupling requires a single rigid bus link and one compliant attachment."))
    link = sc.root
    link.rw_assembly.n_wheels == 0 && isempty(link.faults.faults) &&
        all(name -> isempty(getfield(link, name)), (:thrusters, :magnets, :magnetorquers,
            :laser_terminals, :hinges, :gimbals, :custom_actuators, :solar_arrays,
            :batteries, :power_loads, :sensors, :data_stores, :data_sources)) ||
        throw(ArgumentError("Initial contact coupling requires a bus without mounted actuators, sensors, power/data devices or faults."))
    a = only(sc.attachments)
    a.link === sc.root && length(a.model.bodies) == 1 && length(a.model.joints) == 1 ||
        throw(ArgumentError("Contact coupling requires one body connected directly to the bus mount."))
    a.rest_schedule === nothing && isempty(a.joint_actuators) ||
        throw(ArgumentError("Contact coupling requires fixed joint rest state and no actuators."))
    joint = only(a.model.joints)
    joint.parent == 0 && joint.child == 1 || throw(ArgumentError("Contact coupling requires a mount-to-body joint."))
    for K in (joint.k_translation_n_m, joint.c_translation_n_s_m,
            joint.k_rotation_n_m_rad, joint.c_rotation_n_m_s_rad)
        all(isfinite, K) && K[1,1] >= 0 && K == K[1,1] * SMatrix{3,3,Float64}(I) ||
            throw(ArgumentError("Initial contact coupling requires nonnegative isotropic joint coefficients."))
    end
    CLS.contact_preflight(only(args.contact_models), _contact_bodies(args))
    return nothing
end

function _initialize_contact_runtime!(p, initial_state, owner)
    isempty(p.args.contact_models) && return nothing
    _rhs_env_config(p).execution_mode in (:auto, :serial) ||
        throw(ArgumentError("Contact coupling refuses flat or parallel RHS execution modes."))
    sc = only(initial_state.sc)
    runtime = CLS.contact_prepare(only(p.args.contact_models), _contact_bodies(p.args))
    runtime isa CLS.AbstractContactRuntime || throw(ArgumentError("Contact prepare must return AbstractContactRuntime."))
    run = CLS.ContactRunState(runtime, SVector{3,Float64}(sc.pos), SVector{3,Float64}(sc.vel), Threads.Atomic{Bool}(false))
    owner[] = run
    p.shared_buffers.contact_runtime[] = run
    p.shared_buffers.rhs_plan_override[] = (mode=:serial, allotment=1, scheduler=:serial,
        dominant_axis=:serial, policy_applied=true, effector_decision=_CALIB_SERIAL_EFFECTOR_DECISION)
    return nothing
end

function _apply_contact_loads!(du_view, sc, p, t, forces, torques)
    run = p.shared_buffers.contact_runtime[]
    run === nothing && return nothing
    rt = only(p.shared_buffers.attachment_runtimes)
    CM = SimulationModel.ClothMultibody
    q = SVector{4,Float64}(sc.q); q = q / norm(q)
    aq = SVector{4,Float64}(sc.att_q[:,1]); aq = aq / norm(aq)
    r = SVector{3,Float64}(sc.pos) - run.origin_m - t * run.origin_velocity_m_s
    v = SVector{3,Float64}(sc.vel) - run.origin_velocity_m_s
    states = (CLS.ContactState(r, v, q, SVector{3,Float64}(sc.ω),
            SVector{3,Float64}(forces), CM._rot(q) * SVector{3,Float64}(torques)),
        CLS.ContactState(r + SVector{3,Float64}(sc.att_r[:,1]), v + SVector{3,Float64}(sc.att_v[:,1]),
            aq, SVector{3,Float64}(sc.att_ω[:,1]), SVector{3,Float64}(rt.forces[:,1]),
            CM._rot(aq) * SVector{3,Float64}(rt.torques[:,1])))
    all(s -> all(name -> all(isfinite, getfield(s,name)), fieldnames(CLS.ContactState)), states) && isfinite(t) ||
        throw(ArgumentError("Nonfinite contact trial state or smooth load."))
    loads = CLS.contact_evaluate!(run.runtime, states, Float64(t))
    loads isa NTuple{2,CLS.ContactWrench} || throw(ArgumentError("Contact evaluator must return two COM wrenches."))
    all(w -> all(isfinite,w.force_n) && all(isfinite,w.torque_nm), loads) ||
        throw(ArgumentError("Nonfinite contact wrench."))
    # Closed two-body contact must balance at one common origin. Scale the
    # check by the actual loads, without cancellation at orbital positions.
    scale = max(1.0, norm(loads[1].force_n), norm(loads[2].force_n))
    moment = loads[1].torque_nm + loads[2].torque_nm + cross(SVector{3,Float64}(sc.att_r[:,1]), loads[2].force_n)
    norm(loads[1].force_n + loads[2].force_n) <= 1e-10 * scale &&
        norm(moment) <= 1e-10 * max(scale, norm(loads[1].torque_nm), norm(loads[2].torque_nm)) ||
        throw(ArgumentError("Contact provider violated internal force/moment balance."))
    forces .+= loads[1].force_n
    torques .+= CM._rot(q)' * loads[1].torque_nm
    rt.forces[:,1] .+= loads[2].force_n
    rt.torques[:,1] .+= CM._rot(aq)' * loads[2].torque_nm
    SimulationModel.CompliantAttachmentDynamics.attachment_derivatives!(du_view, sc, rt, 1)
    return nothing
end
