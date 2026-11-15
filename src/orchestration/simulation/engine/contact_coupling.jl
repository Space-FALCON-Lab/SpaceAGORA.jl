# Experimental contact-load orchestration. No provider/native code belongs here.
const CLS = ContactLoadServices

function CLS.with_contact_model(args::SimulationConfiguration, model::CLS.AbstractContactModel)
    return SimulationModel.SimConfig._with_configuration(args; contact_models=(model,))
end

function _contact_bodies(args)
    sc = only(args.dynamics_model.spacecraft)
    bus=CLS.ContactBody(sc.dry_mass + sc.prop_mass, sc.inertia_tensor)
    arm=_robot_arm_coupling(args,1,0.0;run_start=true)
    if arm !== nothing
        CA=SimulationModel.ClothRobotArmDynamics
        return (bus, (CLS.ContactBody(b.mass_kg,CA._link_inertia(b)) for b in arm.plan.model.links)...)
    end
    return (bus, (CLS.ContactBody(b.mass_kg,b.inertia_body_kg_m2)
        for a in sc.attachments for b in a.model.bodies)...)
end

function _validate_contact_coupling!(args, solver_mode)
    isempty(args.contact_models) && return nothing
    length(args.contact_models) == 1 && only(args.contact_models) isa CLS.AbstractContactModel ||
        throw(ArgumentError("Contact coupling requires one AbstractContactModel."))
    length(args.dynamics_model.spacecraft) == 1 || throw(ArgumentError("Contact coupling supports one rigid bus with compliant bodies."))
    solver_mode === :tsit5 && !_run_solver_config(args).parallel ||
        throw(ArgumentError("Contact coupling requires serial :tsit5 integration."))
    isempty(args.external_propagators) && args.interlink_model === nothing ||
        throw(ArgumentError("Contact coupling cannot be combined with external propagation or interlinks."))
    all(e->e isa SimulationModel.RobotArmReactionEffector,args.dynamics_model.dynamic_effectors) && isempty(args.control_model.control_effectors) &&
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
    !_spacecraft_articulated(sc) && length(sc.links) == 1 ||
        throw(ArgumentError("Contact coupling requires a single rigid bus link; articulated backbones are not qualified."))
    link = sc.root
    link.rw_assembly.n_wheels == 0 && isempty(link.faults.faults) &&
        all(name -> isempty(getfield(link, name)), (:thrusters, :magnets, :magnetorquers,
            :laser_terminals, :hinges, :gimbals, :custom_actuators, :solar_arrays,
            :batteries, :power_loads, :sensors, :data_stores, :data_sources)) ||
        throw(ArgumentError("Initial contact coupling requires a bus without mounted actuators, sensors, power/data devices or faults."))
    arm=_robot_arm_coupling(args,1,0.0;run_start=true)
    matrices=Any[]
    if arm !== nothing
        isempty(sc.attachments) && length(args.dynamics_model.dynamic_effectors)==1 ||
            throw(ArgumentError("Contact requires either generic attachments or one compliant robot-arm route, without duplicate arm effectors."))
        !isempty(arm.plan.model.links) && all(arm.plan.q_ref .== arm.plan.q_ref[:,1]) &&
            all(iszero,arm.plan.dq_ref) && all(iszero,arm.plan.ddq_ref) && isempty(arm.joint_actuators) ||
            throw(ArgumentError("Initial arm contact requires a nonempty arm with fixed rest plan and no actuators."))
        CA=SimulationModel.ClothRobotArmDynamics
        for field in (:k_translation_n_m,:c_translation_n_s_m,:k_rotation_n_m_rad,:c_rotation_n_m_s_rad)
            append!(matrices,CA._joint_compliance_matrices(getproperty(arm,field),length(arm.plan.model.links),field))
        end
    else
        isempty(args.dynamics_model.dynamic_effectors) && !isempty(sc.attachments) ||
            throw(ArgumentError("Contact requires compliant attachments or an active compliant robot-arm plan."))
        for a in sc.attachments
            a.link===sc.root && !isempty(a.model.bodies) && a.rest_schedule===nothing && isempty(a.joint_actuators) ||
                throw(ArgumentError("Contact attachments must be mounted on the rigid bus, with fixed rest state and no actuators."))
            for joint in a.model.joints
                append!(matrices,(joint.k_translation_n_m,joint.c_translation_n_s_m,
                    joint.k_rotation_n_m_rad,joint.c_rotation_n_m_s_rad))
            end
        end
    end
    for K in matrices
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

function _apply_contact_loads!(du_view, sc, p, t, forces, torques, arm_loads=nothing)
    run = p.shared_buffers.contact_runtime[]
    run === nothing && return nothing
    rt = only(p.shared_buffers.attachment_runtimes)
    CM = SimulationModel.ClothMultibody
    q = SVector{4,Float64}(sc.q); q = q / norm(q)
    r = SVector{3,Float64}(sc.pos) - run.origin_m - t * run.origin_velocity_m_s
    v = SVector{3,Float64}(sc.vel) - run.origin_velocity_m_s
    bus=CLS.ContactState(r,v,q,SVector{3,Float64}(sc.ω),SVector{3,Float64}(forces),CM._rot(q)*SVector{3,Float64}(torques))
    observations=CLS.ContactState[bus]
    n=arm_loads===nothing ? rt.n_bodies : length(arm_loads.plan.model.links)
    for i in 1:n
        aq=SVector{4,Float64}(arm_loads===nothing ? sc.att_q[:,i] : sc.arm_q[:,i]);aq/=norm(aq)
        if arm_loads===nothing
            push!(observations,CLS.ContactState(r+SVector{3,Float64}(sc.att_r[:,i]),
                v+SVector{3,Float64}(sc.att_v[:,i]),aq,SVector{3,Float64}(sc.att_ω[:,i]),
                SVector{3,Float64}(rt.forces[:,i]),CM._rot(aq)*SVector{3,Float64}(rt.torques[:,i])))
        else
            push!(observations,CLS.ContactState(r+SVector{3,Float64}(sc.arm_r[:,i]-sc.pos),
                SVector{3,Float64}(sc.arm_v[:,i])-run.origin_velocity_m_s,aq,SVector{3,Float64}(sc.arm_ω[:,i]),
                arm_loads.loads.forces[i],CM._rot(aq)*arm_loads.loads.torques_body[i]))
        end
    end
    states=Tuple(observations)
    all(s -> all(name -> all(isfinite, getfield(s,name)), fieldnames(CLS.ContactState)), states) && isfinite(t) ||
        throw(ArgumentError("Nonfinite contact trial state or smooth load."))
    loads = CLS.contact_evaluate!(run.runtime, states, Float64(t))
    loads isa Tuple && length(loads)==length(states) && all(w->w isa CLS.ContactWrench,loads) ||
        throw(ArgumentError("Contact evaluator must return one COM wrench per body in state order."))
    all(w -> all(isfinite,w.force_n) && all(isfinite,w.torque_nm), loads) ||
        throw(ArgumentError("Nonfinite contact wrench."))
    # Closed two-body contact must balance at one common origin. Scale the
    # check by the actual loads, without cancellation at orbital positions.
    scale = max(1.0,sum(w->norm(w.force_n),loads))
    moment=sum(w.torque_nm+cross(s.position_m-r,w.force_n) for (w,s) in zip(loads,states))
    norm(sum(w.force_n for w in loads)) <= 1e-10 * scale &&
        norm(moment) <= 1e-10 * max(scale,sum(w->norm(w.torque_nm),loads)) ||
        throw(ArgumentError("Contact provider violated internal force/moment balance."))
    forces .+= loads[1].force_n
    torques .+= CM._rot(q)' * loads[1].torque_nm
    for i in 1:n
        w=loads[i+1]; tau=CM._rot(states[i+1].quaternion)'*w.torque_nm
        if arm_loads===nothing
            rt.forces[:,i] .+= w.force_n
            rt.torques[:,i] .+= tau
        else
            arm_loads.loads.forces[i]+=w.force_n
            arm_loads.loads.torques_body[i]+=tau
        end
    end
    if arm_loads===nothing
        for k in 1:rt.n_att
            SimulationModel.CompliantAttachmentDynamics.attachment_derivatives!(du_view,sc,rt,k)
        end
    else
        SimulationModel.ClothRobotArmDynamics._coupled_robot_arm_derivatives!(du_view,sc,arm_loads.plan,arm_loads.loads)
    end
    return nothing
end
