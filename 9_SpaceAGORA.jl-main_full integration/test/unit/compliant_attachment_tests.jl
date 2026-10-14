module CompliantAttachmentTests

# Compliant multibody attachments mounted on spacecraft links (CompliantAttachment): state,
# loads, two-way coupling through the rigid bus and the articulated backbone.

using Test
using LinearAlgebra
using StaticArrays
using DataFrames
using SpaceAGORA
using SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const AB = SM.ArticulatedBody
const CM = SM.ClothMultibody
const CAD = SM.CompliantAttachmentDynamics
const SE = SpaceAGORA.SimulationEngine

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

rotq(axis, ang) = (a = normalize(collect(Float64, axis)); SVector{4, Float64}(a[1] * sin(ang / 2), a[2] * sin(ang / 2), a[3] * sin(ang / 2), cos(ang / 2)))
Rmat(q) = AB._rotmat(SVector{4, Float64}(q))
const QI = SVector(0.0, 0.0, 0.0, 1.0)
const Z3 = SVector(0.0, 0.0, 0.0)

function mklink(; root=false, m, dims, r=(0.0, 0.0, 0.0), q=(0.0, 0.0, 0.0, 1.0), kwargs...)
    return SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(r...), q=MVector{4, Float64}(q...); kwargs...)
end

function mkjoint(par::SM.Link, child::SM.Link, P; kwargs...)
    rp = par.root ? zeros(3) : collect(par.r)
    Rp = par.root ? Matrix(1.0I, 3, 3) : Matrix(Rmat(par.q))
    Rc = Matrix(Rmat(child.q))
    p1 = SVector{3, Float64}(Rp' * (collect(P) - rp))
    p2 = SVector{3, Float64}(Rc' * (collect(P) - collect(child.r)))
    return SM.Joint(par, p1, child, p2; kwargs...)
end

const T0 = SM.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)

# Allocation counts are only meaningful without coverage instrumentation, which
# blocks inlining and adds allocations (cf. test/unit/parallel/cost_robust_timing_tests.jl).
const _ALLOC_CHECKS = Base.JLOptions().code_coverage == 0

# The attachment state holds ABSOLUTE inertial positions (~7e6 m here, ulp 1e-9 m), so the spring forces carry a
# roundoff floor of k * 1e-9 m. Absolute tolerances below that floor only shrink the step; `abs_orbit` stays at it.
function tight_tolerances(; dt_max=0.05, rel=1e-12, abs_orbit=1e-9, abs_att=1e-9)
    return SM.IntegrationTolerances(
        reltol_orbit=rel, abstol_orbit=abs_orbit, reltol_atmosphere=rel, abstol_atmosphere=abs_orbit,
        reltol_quaternion=rel, abstol_quaternion=abs_att, reltol_mass=rel, abstol_mass=abs_att,
        reltol_angular_rate=rel, abstol_angular_rate=abs_att, dt_max_orbit=dt_max, dt_max_atmosphere=dt_max,
    )
end

function engine_args(sc; mission_time, data_rate=0.05, effectors=(), control=(), control_rates=Float64[],
        orientation=true, solver_mode=:dp8, tol=tight_tolerances(), results_directory=mktempdir())
    scs = sc isa AbstractVector ? sc : [sc]
    base = make_example_config(
        planet=SM.make_no_gram_planet(:earth), spacecraft=scs[1], mission_time=mission_time, initial_time=T0,
        dynamic_effectors=effectors, density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=orientation, keplerian=true, verbose=false, results=false, results_directory=results_directory,
        solver_config=SM.SolverConfig(solver_mode=solver_mode),
    )
    return SM.SimConfig._with_configuration(base;
        integration_tolerances=tol,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=orientation, num_steps_to_save=1000, data_rate=data_rate),
        control_model=SM.ControlModel(control_effectors=control, control_rates=control_rates),
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[scs...], effectors),
    )
end

ic_at(r, v; q=QI, ω=Z3) = SM.CartesianInitialCondition(SVector{3, Float64}(r), SVector{3, Float64}(v); q=SVector{4, Float64}(q), ang_vel=SVector{3, Float64}(ω))

# Small compliant chain along the mount-frame x axis: `n` bodies of mass `m`, joint 1 attaches body 1 to the
# mount, joint j attaches body j to body j - 1. Isotropic springs (kx N/m, kr N m/rad) and dampers (cx, cr).
# The initial state is deflected (displaced, rotated, moving bodies); every rest orientation is the identity.
function chain_build(; n=3, m=1.0, kx=200.0, cx=0.0, kr=20.0, cr=0.0, deflect=true)
    J = SM.thin_panel_inertia(m, 1.0, 0.4)
    nodes = SM.CompliantTopologyNode[]
    for i in 1:n
        pos = SVector(i - 0.5, 0.0, 0.0)
        q = QI; v = Z3; w = Z3
        if deflect && i == n
            pos += SVector(0.0, 0.05, 0.02); v = SVector(0.0, 0.1, 0.0)
        end
        if deflect && i == 2
            q = rotq([0, 0, 1], 0.1); w = SVector(0.2, 0.0, 0.1)
        end
        push!(nodes, SM.CompliantTopologyNode(Symbol("b$i"); mass_kg=m, inertia_body_kg_m2=J, position=pos, quaternion=q, velocity=v, angular_velocity=w))
    end
    edges = SM.CompliantTopologyEdge[]
    for i in 1:n
        push!(edges, SM.CompliantTopologyEdge(Symbol("j$i"), i - 1, i;
            parent_point_body=i == 1 ? Z3 : SVector(0.5, 0.0, 0.0), child_point_body=SVector(-0.5, 0.0, 0.0),
            k_translation_n_m=kx, c_translation_n_s_m=cx, k_rotation_n_m_rad=kr, c_rotation_n_m_s_rad=cr,
            rest_child_parent_quat=QI))
    end
    return SM.build_compliant_topology(nodes, edges)
end

# Bus-like rigid spacecraft carrying attachments.
function rigid_sc(atts; m=20.0, dims=(1.0, 0.8, 0.6), ic=ic_at((7.0e6, 0.0, 0.0), (0.0, 0.0, 0.0)), bus=nothing)
    bus === nothing && (bus = mklink(root=true, m=m, dims=dims))
    return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=ic, inertia_tensor=bus.inertia,
        attachments=atts isa Function ? atts(bus) : atts), bus
end

# bus + one hinged panel; attachments are built from the panel link
function hinge_sc(make_atts; k=4.0, c=0.0, θ0=0.3, ic=ic_at((7.0e6, 0.0, 0.0), (30.0, -20.0, 10.0)))
    bus = mklink(root=true, m=10.0, dims=(1.0, 1.0, 1.0))
    panel = mklink(m=2.0, dims=(0.05, 1.0, 0.5), r=(0.0, 1.1, 0.0))
    joint = mkjoint(bus, panel, (0.0, 0.5, 0.0); joint_type=:hinge, axis=[0, 0, 1], stiffness=k, damping=c, initial_q=θ0)
    return SM.SpacecraftModel(; joints=[joint], links=[bus, panel], root=bus, initial_condition=ic, attachments=make_atts(bus, panel)), bus, panel
end

function column(df, name, k)
    return df[!, "sc1_$(name)_$k"]
end

# ---------------------------------------------------------------------------
# (a) no attachments: nothing changes
# ---------------------------------------------------------------------------

@testset "no attachments: runs are bit-identical and runtimes are empty" begin
    function build(kw)
        bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
        L = mklink(m=2.0, dims=(0.05, 1.0, 0.5), r=(0.0, -1.1, 0.0))
        return SM.SpacecraftModel(; joints=[mkjoint(bus, L, (0.0, -0.5, 0.0))], links=[bus, L], root=bus,
            initial_condition=ic_at((7.0e6, 0.0, 0.0), (0.0, 7.5e3, 0.0); ω=(0.0, 0.0, 1e-3)), inertia_tensor=bus.inertia, kw...)
    end
    run_it(sc) = SpaceAGORA.run_simulation(engine_args(sc; mission_time=60.0, data_rate=5.0, solver_mode=:tsit5,
        effectors=(SM.InverseSquaredJ2GravityModel(),), tol=SM.IntegrationTolerances()); return_results=true)
    omitted = run_it(build(()))
    empty = run_it(build((attachments=SM.CompliantAttachment[],)))
    @test names(omitted.table) == names(empty.table)
    @test all(isequal(omitted.table[!, n], empty.table[!, n]) for n in names(omitted.table))
    @test !any(n -> occursin("attachment", n), names(empty.table))
    @test !("sc1_system_com_1" in names(empty.table))
    sc = build((attachments=SM.CompliantAttachment[],))
    @test isempty(sc.attachments)
    args = engine_args(sc; mission_time=1.0)
    p = SM.ODEParams(n_sats=1, args=args)
    SE._initialize_attachment_runtimes!(p)
    @test all(isnothing, p.shared_buffers.attachment_runtimes)
    @test !p.shared_buffers.attachments_present[]
    u = SE.build_initial_conditions(args)
    @test !hasproperty(u.sc[1], :att_r)

    # An articulated spacecraft with attachments=[] equals the one without the keyword.
    function art(kw)
        bus = mklink(root=true, m=10.0, dims=(1.0, 1.0, 1.0))
        panel = mklink(m=2.0, dims=(0.05, 1.0, 0.5), r=(0.0, 1.1, 0.0))
        joint = mkjoint(bus, panel, (0.0, 0.5, 0.0); joint_type=:hinge, axis=[0, 0, 1], stiffness=4.0, damping=0.1, initial_q=0.3)
        return SM.SpacecraftModel(; joints=[joint], links=[bus, panel], root=bus, initial_condition=ic_at((7.0e6, 0.0, 0.0), (30.0, -20.0, 10.0)), kw...)
    end
    a1 = SpaceAGORA.run_simulation(engine_args(art(()); mission_time=20.0, data_rate=1.0, effectors=(SM.InverseSquaredGravityModel(),)); return_results=true)
    a2 = SpaceAGORA.run_simulation(engine_args(art((attachments=SM.CompliantAttachment[],)); mission_time=20.0, data_rate=1.0, effectors=(SM.InverseSquaredGravityModel(),)); return_results=true)
    @test names(a1.table) == names(a2.table)
    @test all(isequal(a1.table[!, n], a2.table[!, n]) for n in names(a1.table))
end

# ---------------------------------------------------------------------------
# (b) standalone equivalence
# ---------------------------------------------------------------------------

# RK4 on the standalone model with the exact time-dependent rest quaternions at each stage.
function standalone_rk4(model, x0, rest_at, actuators, dt, nsteps)
    f(x, t) = CM.compliant_multibody_dynamics(model, x, t; joint_rest_quaternions=rest_at(t), joint_actuators=actuators)
    xs = Vector{Vector{Float64}}(undef, nsteps + 1)
    xs[1] = copy(x0)
    x = copy(x0)
    for k in 1:nsteps
        t = (k - 1) * dt
        k1 = f(x, t); k2 = f(x + 0.5dt * k1, t + 0.5dt); k3 = f(x + 0.5dt * k2, t + 0.5dt); k4 = f(x + dt * k3, t + dt)
        x = x + dt / 6 * (k1 + 2k2 + 2k3 + k4)
        CM._normalize_state_quaternions!(x)
        xs[k + 1] = copy(x)
    end
    return xs
end

function max_relative_position_error(table, xs, ratio, nb, mount_r, R_m)
    # engine positions are inertial; the standalone ones are in the mount frame of a fixed mount
    err = 0.0
    for row in 1:nrow(table)
        k = round(Int, table.time[row] / 0.05) * ratio + 1
        for i in 1:nb
            eng = [table[row, "sc1_attachment_pose_$(7i - 7 + c)"] for c in 1:3]
            ref = mount_r + R_m * collect(xs[k][(13 * (i - 1) + 1):(13 * (i - 1) + 3)])
            err = max(err, norm(eng - ref))
        end
    end
    return err
end

@testset "engine attachment matches the standalone model on a (nearly) fixed mount" begin
    build = chain_build(n=3, kx=200.0, cx=2.0, kr=20.0, cr=0.5)
    mp = SVector(0.5, 0.2, 0.1)
    mq = rotq([0, 0, 1], 0.4)
    bus_q = rotq([1, 2, 3], 0.3)
    tend = 3.0
    dt = 1.0e-3
    mass_bus = 1.0e8
    for (label, schedule) in (("constant rest", nothing), ("time-varying rest", (out, t) -> (out[2] = SVector(0.0, 0.0, sin(0.5 * 0.1 * sin(2t)), cos(0.5 * 0.1 * sin(2t))); nothing)))
        act = [CM.CompliantJointActuator(:a2, 2; kp_n_m_rad=3.0, kd_n_m_s_rad=0.2, feedforward_torque_child_body=(0.0, 0.0, 0.05), torque_limit_n_m=1.0)]
        function make_sc()
            bus = mklink(root=true, m=mass_bus, dims=(1.0, 1.0, 1.0))
            att = SM.CompliantAttachment(; model=build, link=bus, mount_point=mp, mount_quaternion=mq, joint_actuators=act, rest_schedule=schedule)
            return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=ic_at((7.0e6, 0.0, 0.0), (0.0, 0.0, 0.0); q=bus_q),
                inertia_tensor=bus.inertia, attachments=[att])
        end
        res = SpaceAGORA.run_simulation(engine_args(make_sc(); mission_time=tend, data_rate=0.05, tol=tight_tolerances(dt_max=0.02, rel=1e-12)); return_results=true)
        rest0 = [j.rest_child_parent_quat for j in build.model.joints]
        rest_at = t -> begin
            out = copy(rest0)
            schedule === nothing || schedule(out, t)
            out
        end
        xs = standalone_rk4(CM.CompliantMultibodyModel(build.model.bodies, build.model.joints, Z3, QI), build.initial_state,
            rest_at, act, dt, round(Int, tend / dt))
        mount_r = SVector(7.0e6, 0.0, 0.0) + Rmat(bus_q) * mp
        R_m = Rmat(SVector{4, Float64}(CM._quat_mul(bus_q, mq)))
        err = max_relative_position_error(res.table, xs, 50, 3, mount_r, R_m)
        # the bus is the only physical difference: it recoils by about m_att/m_bus x deflection (~1e-9 m); the integrator
        # tolerance floor is reltol*|r| = 1e-12 * 7e6 = 7e-6 m, and the standalone RK4 error is ~1e-12 m. Measured 2.4e-7 m.
        @info "standalone vs engine ($label)" max_position_error_m = err
        @test err < 2e-6
        # attachment_pose is inertial (pos + att_r)
    # attachment_pose has 7 numbers per body and unit quaternions
        @test count(startswith("sc1_attachment_pose_"), names(res.table)) == 21
        @test all(abs(norm([res.table[end, "sc1_attachment_pose_$(7i - 7 + c)"] for c in 4:7]) - 1) < 1e-9 for i in 1:3)
        # the bus barely moves
        @test norm([res.table[end, "sc1_pos_$c"] for c in 1:3] - [7.0e6, 0, 0]) < 1e-5
    end
end

# ---------------------------------------------------------------------------
# (c) two-way coupling: conservation
# ---------------------------------------------------------------------------

# Bodies of spacecraft + attachments: (mass, position, velocity, attitude quaternion, body rate, body-frame inertia).
function system_bodies(sc, tree, sv)
    bodies = Tuple{Float64, SVector{3, Float64}, SVector{3, Float64}, SVector{4, Float64}, SVector{3, Float64}, SMatrix{3, 3, Float64, 9}}[]
    kin = nothing
    if tree === nothing
        push!(bodies, (sv.mass, SVector{3}(sv.pos), SVector{3}(sv.vel), SVector{4}(sv.q), SVector{3}(sv.ω), SMatrix{3, 3, Float64, 9}(sc.inertia_tensor)))
    else
        base = AB.ArticulatedBaseState(sv.pos, sv.vel, sv.q, sv.ω)
        kin = AB.articulated_kinematics(tree, base, sv.joint_q, sv.joint_qd)
        masses = copy(tree.mass); masses[1] = sv.mass - AB.articulated_moving_mass(tree)
        for b in 1:tree.nb
            push!(bodies, (masses[b], kin.pos[b], kin.vel[b], kin.quat[b], kin.ω[b], tree.inertia[b]))
        end
    end
    rt = CAD.build_attachment_runtime(sc, tree)
    pe = 0.0
    if rt !== nothing
        for k in 1:rt.n_att
            a = rt.attachments[k]
            b = rt.link_body[k]
            mount = if tree === nothing
                CM.compliant_mount_kinematics(sv.pos, sv.vel, sv.q, sv.ω, rt.link_offset[k], rt.link_q[k], a.mount_point, a.mount_quaternion)
            else
                CM.compliant_mount_kinematics(kin.pos[b], kin.vel[b], kin.quat[b], kin.ω[b], rt.link_offset[k], rt.link_q[k], a.mount_point, a.mount_quaternion)
            end
            col(i) = rt.col0[k] + i
            # att_r/att_v are relative to the spacecraft pos/vel; the invariants use inertial quantities
            st(i) = (r=SVector{3}(sv.pos) + SVector{3}(sv.att_r[:, col(i)]), q=SVector{4}(sv.att_q[:, col(i)]), v=SVector{3}(sv.vel) + SVector{3}(sv.att_v[:, col(i)]), ω=SVector{3}(sv.att_ω[:, col(i)]))
            for (i, body) in pairs(a.model.bodies)
                s = st(i)
                push!(bodies, (body.mass_kg, s.r, s.v, s.q, s.ω, SMatrix{3, 3, Float64, 9}(body.inertia_body_kg_m2)))
            end
            # compliant spring energy: 1/2 dr' K dr + 1/2 phi' K phi, exact for isotropic stiffness
            for (jidx, j) in pairs(a.model.joints)
                pq = j.parent == 0 ? mount.q : st(j.parent).q
                pr = j.parent == 0 ? mount.r : st(j.parent).r
                ppoint = pr + CM._rot(pq) * j.parent_point_body
                cs = st(j.child)
                cpoint = cs.r + CM._rot(cs.q) * j.child_point_body
                Δr = ppoint - cpoint
                qerr = CM._quat_mul(CM._quat_mul(pq, rt.rest[k][jidx]), CM._quat_conj(cs.q))
                ϕ = CM._axis_angle_error(qerr)
                pe += dot(Δr, j.k_translation_n_m * Δr) / 2 + dot(ϕ, j.k_rotation_n_m_rad * ϕ) / 2
            end
        end
    end
    tree === nothing || (pe += AB.articulated_potential_energy(tree, sv.joint_q))
    return bodies, pe
end

# Positions enter relative to the first body, so the measurement itself carries no 7e6 m cancellation.
function invariants(sc, tree, sv)
    bodies, pe = system_bodies(sc, tree, sv)
    M = sum(b[1] for b in bodies)
    r0 = bodies[1][2]
    dX = sum(b[1] * (b[2] - r0) for b in bodies) / M
    V = sum(b[1] * b[3] for b in bodies) / M
    L = zero(SVector{3, Float64}); KE = 0.0
    for (m, r, v, q, ω, J) in bodies
        L += m * cross(r - r0 - dX, v - V) + Rmat(q) * (J * ω)
        KE += m * dot(v, v) / 2 + dot(ω, J * ω) / 2
    end
    return M * V, L, KE + pe
end

function conservation_drifts(sc; mission_time, solver_mode=:dp8)
    tree = SM.articulated_has_moving_joints(sc) ? AB.build_articulated_tree(sc) : nothing
    args = engine_args(sc; mission_time=mission_time, data_rate=0.1, solver_mode=solver_mode, tol=tight_tolerances(dt_max=0.02))
    res = SpaceAGORA.run_simulation(args; return_results=true, return_solution=true)
    sol = res.solution
    P0, L0, E0 = invariants(sc, tree, sol.u[1].sc[1])
    dP = 0.0; dL = 0.0; dE = 0.0
    for u in sol.u
        P, L, E = invariants(sc, tree, u.sc[1])
        dP = max(dP, norm(P - P0) / norm(P0)); dL = max(dL, norm(L - L0) / norm(L0)); dE = max(dE, abs(E - E0) / abs(E0))
    end
    return (dP=dP, dL=dL, dE=dE, sol=sol, res=res, steps=length(sol.t), E0=E0)
end

@testset "free-floating bus + attachment conserves momentum, angular momentum and energy" begin
    build = chain_build(n=3, kx=200.0, cx=0.0, kr=20.0, cr=0.0)
    ic = ic_at((7.0e6, 0.0, 0.0), (30.0, -20.0, 10.0); q=rotq([1, -1, 2], 0.9), ω=(0.03, -0.02, 0.05))
    sc, _ = rigid_sc(bus -> [SM.CompliantAttachment(; model=build, link=bus, mount_point=(0.6, 0.1, -0.2), mount_quaternion=rotq([1, 1, 0], 0.5))]; ic=ic)
    d = conservation_drifts(sc; mission_time=20.0)
    @info "rigid bus + attachment conservation (20 s)" linear_momentum = d.dP angular_momentum = d.dL energy = d.dE steps = d.steps
    @test d.dP < 1e-13
    # measured 2e-10 (angular momentum), 1.4e-12 (energy); set by the 1e-9 absolute tolerances over |L0| ~ 0.1
    @test d.dL < 1e-8
    @test d.dE < 1e-10
    # the attachment really moves relative to the bus: relative position of the last body to the bus changes by O(deflection)
    r_last = [d.res.table[:, "sc1_attachment_pose_$(14 + c)"] .- d.res.table[:, "sc1_pos_$c"] for c in 1:3]
    @test maximum(abs, r_last[2] .- r_last[2][1]) > 1e-2
    # bus mass is the spacecraft mass only; the system total is the spacecraft plus the attachments
    @test d.res.table.sc1_mass[end] == 20.0
    @test SM.attachment_total_mass(sc) == 3.0

    # attachment on a moving hinged link of an articulated spacecraft: the reaction goes through the backbone wrench path
    sca, _, _ = hinge_sc((bus, panel) -> [SM.CompliantAttachment(; model=build, link=panel, mount_point=(0.0, 0.5, 0.0), mount_quaternion=rotq([1, 0, 0], 0.7))]; k=4.0, c=0.0, θ0=0.3)
    da = conservation_drifts(sca; mission_time=20.0)
    @info "articulated (hinged link) + attachment conservation (20 s)" linear_momentum = da.dP angular_momentum = da.dL energy = da.dE steps = da.steps
    @test da.dP < 1e-13
    @test da.dL < 1e-8                                                  # measured 2e-9
    @test da.dE < 1e-10
    @test abs(da.sol.u[end].sc[1].joint_q[1] - 0.3) > 1e-3             # the hinge moved
    # system_com includes the attachment bodies: it moves at the (conserved) initial COM velocity
    tbl = da.res.table
    P0, _, _ = invariants(sca, AB.build_articulated_tree(sca), da.sol.u[1].sc[1])
    v0 = collect(P0) / (da.sol.u[1].sc[1].mass + SM.attachment_total_mass(sca))     # initial system COM velocity
    com_rate = [(tbl[end, "sc1_system_com_$i"] - tbl[1, "sc1_system_com_$i"]) / (tbl.time[end] - tbl.time[1]) for i in 1:3]
    @test norm(com_rate - v0) < 1e-3
    # the same invariant check but with the attachment on the root body (articulated spacecraft, root mount)
    scr, _, _ = hinge_sc((bus, panel) -> [SM.CompliantAttachment(; model=build, link=bus, mount_point=(0.6, 0.0, 0.0))]; k=4.0, c=0.0, θ0=0.3)
    dr = conservation_drifts(scr; mission_time=10.0)
    @info "articulated (root mount) + attachment conservation (10 s)" linear_momentum = dr.dP angular_momentum = dr.dL energy = dr.dE
    @test dr.dP < 1e-13 && dr.dL < 1e-8 && dr.dE < 1e-10
end

# ---------------------------------------------------------------------------
# (d) rigid limit: stiff damped attachment body vs. the same body as a :fixed link
# ---------------------------------------------------------------------------

# The attachment state is relative to the bus, so spring forces carry no k * ulp(7e6 m) roundoff floor and a stiff
# attachment (2e5 N/m, 316 rad/s) integrates with ordinary tolerances (this case failed with MaxIters when it was absolute).
@testset "rigid limit: stiff attachments vs fixed links, circular orbit" begin
    planet = SM.make_no_gram_planet(:earth)
    r0 = 7.0e6
    vc = sqrt(planet.μ / r0)
    n = vc / r0
    ic = ic_at((r0, 0.0, 0.0), (0.0, vc, 0.0); ω=(0.0, 0.0, n))
    panel_dims = (0.05, 1.0, 0.5)
    # two single-body attachments on the +/- y faces of the bus (symmetric: the system COM is the bus origin)
    function fixed_case()
        bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
        L = mklink(m=2.0, dims=panel_dims, r=(0.0, -1.1, 0.0))
        R = mklink(m=2.0, dims=panel_dims, r=(0.0, 1.1, 0.0))
        sc = SM.SpacecraftModel(; joints=[mkjoint(bus, L, (0.0, -0.5, 0.0)), mkjoint(bus, R, (0.0, 0.5, 0.0))], links=[bus, L, R], root=bus, initial_condition=ic)
        sc.inertia_tensor = AB.build_articulated_tree(sc).inertia[1]
        return sc
    end
    function attached_case()
        bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
        panel = mklink(m=2.0, dims=panel_dims)               # only used for its mass and inertia
        function att(sgn)
            kx = 2.0e5; kr = 2.0e4
            cx = 2 * 0.7 * sqrt(kx * 2.0); cr = 2 * 0.7 * sqrt(kr * panel.inertia[1, 1])
            node = SM.CompliantTopologyNode(:panel; mass_kg=2.0, inertia_body_kg_m2=panel.inertia, position=SVector(0.0, sgn * 0.6, 0.0))
            edge = SM.CompliantTopologyEdge(:hold, 0, 1; parent_point_body=Z3, child_point_body=SVector(0.0, -sgn * 0.6, 0.0),
                k_translation_n_m=kx, c_translation_n_s_m=cx, k_rotation_n_m_rad=kr, c_rotation_n_m_s_rad=cr)
            return SM.CompliantAttachment(; model=SM.build_compliant_topology([node], [edge]), link=bus, mount_point=(0.0, sgn * 0.5, 0.0))
        end
        return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=ic, inertia_tensor=SMatrix{3, 3, Float64, 9}(AB.build_articulated_tree(fixed_case()).inertia[1]),
            attachments=[att(-1.0), att(1.0)])
    end
    fixed_sc = fixed_case()
    att_sc = attached_case()
    # mass bookkeeping: the fixed links are in dry_mass, the attachment bodies are not
    @test fixed_sc.dry_mass == 24.0
    @test att_sc.dry_mass == 20.0
    @test SM.attachment_total_mass(att_sc) == 4.0
    tend = 100.0
    run_it(sc, mode) = SpaceAGORA.run_simulation(engine_args(sc; mission_time=tend, data_rate=2.0, effectors=(SM.InverseSquaredGravityModel(),),
        solver_mode=mode, tol=tight_tolerances(dt_max=1.0, rel=1e-12, abs_orbit=1e-9, abs_att=1e-9)); return_results=true)
    rig = run_it(fixed_sc, :dp8)
    t_att = @elapsed att = run_it(att_sc, :dp8)
    @info "stiff attachment run (2e5 N/m), dp8" wall_s = t_att
    @test nrow(att.table) == nrow(rig.table)
    dcom = maximum(norm([att.table[i, "sc1_system_com_$k"] - rig.table[i, "sc1_pos_$k"] for k in 1:3]) for i in 1:nrow(att.table))
    dq = maximum(norm([att.table[i, "sc1_q_$k"] - rig.table[i, "sc1_q_$k"] for k in 1:4]) for i in 1:nrow(att.table))
    @info "rigid limit, stiff attachments vs fixed links (100 s)" system_com_max_diff_m = dcom bus_quaternion_max_diff = dq
    @test dcom < 1e-6            # measured 5.5e-8 m; dominated by the 1e-12 relative tolerance on a 7e6 m position (7e-6 m) and abs_orbit = 1e-7 m
    @test dq < 1e-11             # measured 5.9e-14
    @test att.table.sc1_mass[end] == 20.0 && rig.table.sc1_mass[end] == 24.0
    # the stiff-limit stays rigid: body separation from the bus equals the configured 1.1 m
    sep = [norm([att.table[i, "sc1_attachment_pose_$c"] - att.table[i, "sc1_pos_$c"] for c in 1:3]) for i in 1:nrow(att.table)]
    @test maximum(abs.(sep .- 1.1)) < 1e-5
    # a stiff model also runs under the implicit solver (cloth meshes are stiff in general)
    ras = SpaceAGORA.run_simulation(engine_args(att_sc; mission_time=40.0, data_rate=2.0, effectors=(SM.InverseSquaredGravityModel(),),
        solver_mode=:rodas5p, tol=tight_tolerances(dt_max=2.0, rel=1e-9, abs_orbit=1e-6, abs_att=1e-7)); return_results=true)
    @test ras.table.time[end] ≈ 40.0
    j = findlast(<=(40.0), rig.table.time)
    @test norm([ras.table[end, "sc1_system_com_$k"] - rig.table[j, "sc1_pos_$k"] for k in 1:3]) < 1e-1
end

# ---------------------------------------------------------------------------
# (e) per-body external wrench in the articulated library
# ---------------------------------------------------------------------------

@testset "articulated_dynamics!: per-body wrenches" begin
    # Locked chain: hinge + slide + ball at their rest coordinates, nothing moving. A force m_b a on every
    # body's COM (a force through the system COM overall) must give pure translation: the root accelerates
    # at a, nothing rotates and no joint accelerates.
    bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
    A = mklink(m=3.0, dims=(0.1, 1.2, 0.6), r=(0.0, 1.0, 0.2), q=rotq([1, 0, 0], 0.3))
    B = mklink(m=2.0, dims=(0.2, 0.9, 0.4), r=(0.3, 2.2, 0.3), q=rotq([0, 1, 0], -0.5))
    C = mklink(m=1.5, dims=(0.3, 0.3, 0.6), r=(0.4, 3.0, 0.1), q=rotq([1, 1, 0], 0.8))
    joints = [
        mkjoint(bus, A, (0.0, 0.6, 0.1); joint_type=:hinge, axis=[0.3, 0.2, 0.9], stiffness=6.0),
        mkjoint(A, B, (0.1, 1.7, 0.25); joint_type=:slide, axis=[1.0, 0.5, 0.0], stiffness=15.0),
        mkjoint(B, C, (0.35, 2.7, 0.2); joint_type=:ball, stiffness=3.0, rest=rotq([0, 0, 1], 0.2), initial_q=rotq([0, 0, 1], 0.2)),
    ]
    sc = SM.SpacecraftModel(; joints=joints, links=[bus, A, B, C], root=bus, initial_condition=ic_at((7.0e6, 0.0, 0.0), (0.0, 0.0, 0.0)))
    tree = AB.build_articulated_tree(sc)
    ws = AB.ArticulatedWorkspace(tree)
    jq = copy(tree.q0); jqd = zeros(tree.nv)
    base = AB.ArticulatedBaseState(zeros(3), zeros(3), rotq([1, -2, 1], 0.7), zeros(3))
    nog = r -> Z3
    a_target = SVector(0.7, -0.3, 1.1)
    Fs = [tree.mass[b] * a_target for b in 1:tree.nb]
    Ts = fill(Z3, tree.nb)
    a, α, qdd = AB.articulated_dynamics!(ws, tree, base, jq, jqd, Z3, Z3, nog; body_force_world=Fs, body_torque_world=Ts)
    @info "locked chain, force through every COM" accel_error = norm(a - a_target) angular = norm(α) joint = maximum(abs, qdd)
    @test norm(a - a_target) < 1e-12
    @test norm(α) < 1e-12
    @test maximum(abs, qdd) < 1e-12
    # Without the wrench keywords the call is the Phase-1 call.
    a0, α0, q0 = AB.articulated_dynamics!(ws, tree, base, jq, jqd, Z3, Z3, nog)
    @test norm(a0) < 1e-12 && norm(α0) < 1e-12 && maximum(abs, q0) < 1e-12

    # Pure torque on the panel of a bus + hinge system: [M11 M12; M12 M22] [psi''; theta''] = tau [1; 1].
    hb = mklink(root=true, m=10.0, dims=(1.0, 1.0, 1.0))
    hp = mklink(m=2.0, dims=(0.05, 1.0, 0.5), r=(0.0, 1.1, 0.0))
    hsc = SM.SpacecraftModel(; joints=[mkjoint(hb, hp, (0.0, 0.5, 0.0); joint_type=:hinge, axis=[0, 0, 1], stiffness=0.0)], links=[hb, hp], root=hb,
        initial_condition=ic_at((7.0e6, 0.0, 0.0), (0.0, 0.0, 0.0)))
    ht = AB.build_articulated_tree(hsc)
    hws = AB.ArticulatedWorkspace(ht)
    Ib = 10.0 / 12 * 2.0; Ip = 2.0 / 12 * (0.05^2 + 1.0^2); aa, bb = 0.5, 0.6; μ = 10.0 * 2.0 / 12.0
    M11 = Ib + Ip + μ * (aa + bb)^2; M12 = Ip + μ * bb * (aa + bb); M22 = Ip + μ * bb^2
    τ = 0.37
    expected = [M11 M12; M12 M22] \ [τ, τ]
    hbase = AB.ArticulatedBaseState(zeros(3), zeros(3), QI, zeros(3))
    ha, hα, hq = AB.articulated_dynamics!(hws, ht, hbase, [0.0], [0.0], Z3, Z3, nog;
        body_force_world=[Z3, Z3], body_torque_world=[Z3, SVector(0.0, 0.0, τ)])
    @info "torque on the panel" psi_ddot = hα[3] expected_psi = expected[1] theta_ddot = hq[1] expected_theta = expected[2]
    @test hα[3] ≈ expected[1] rtol = 1e-12
    @test hq[1] ≈ expected[2] rtol = 1e-12
    # A force on the root through its COM equals the base_force_world input.
    F = SVector(1.0, -2.0, 0.5)
    r1 = AB.articulated_dynamics!(hws, ht, hbase, [0.0], [0.0], F, Z3, nog)
    r1 = (collect(r1[1]), collect(r1[2]), copy(r1[3]))
    r2 = AB.articulated_dynamics!(hws, ht, hbase, [0.0], [0.0], Z3, Z3, nog; body_force_world=[F, Z3], body_torque_world=[Z3, Z3])
    @test collect(r2[1]) ≈ r1[1] rtol = 1e-13
    @test collect(r2[2]) ≈ r1[2] atol = 1e-13
    @test r2[3] ≈ r1[3] rtol = 1e-13
    # kinematics in place equals the allocating version
    kin = AB.articulated_kinematics(tree, base, jq, jqd)
    buf = (pos=ws.kpos, quat=ws.kquat, vel=ws.kvel, ω=ws.kwb, ω_world=ws.kww)
    AB.articulated_kinematics!(buf, tree, base, jq, jqd)
    @test buf.pos == kin.pos && buf.quat == kin.quat && buf.vel == kin.vel && buf.ω == kin.ω
    # in-place kinematics are allocation-free
    AB.articulated_kinematics!(buf, tree, base, jq, jqd)
    @test (@allocated AB.articulated_kinematics!(buf, tree, base, jq, jqd)) == 0 skip=!_ALLOC_CHECKS
end

@testset "compliant_joint_loads_in_place! agrees with compliant_joint_loads" begin
    build = chain_build(n=3, kx=200.0, cx=2.0, kr=20.0, cr=0.5)
    model = CM.CompliantMultibodyModel(build.model.bodies, build.model.joints, Z3, QI)
    x = build.initial_state
    rest = [rotq([0, 0, 1], 0.05), rotq([0, 1, 0], -0.02), QI]
    act = [CM.CompliantJointActuator(:a, 2; kp_n_m_rad=3.0, kd_n_m_s_rad=0.2, feedforward_torque_child_body=(0.0, 0.0, 0.05), torque_limit_n_m=1.0)]
    loads = CM.compliant_joint_loads(model, x; joint_rest_quaternions=rest, joint_actuators=act)
    pos = zeros(3, 3); quat = zeros(4, 3); vel = zeros(3, 3); ω = zeros(3, 3)
    for i in 1:3
        s = CM.compliant_state_parts(x, i)
        pos[:, i] .= s.r; quat[:, i] .= s.q; vel[:, i] .= s.v; ω[:, i] .= s.ω
    end
    F = zeros(3, 3); T = zeros(3, 3)
    mount = CM.CompliantMountKinematics(Z3, QI, Z3, Z3)
    restv = SVector{4, Float64}.(rest)
    cmodel = CM.compile_compliant_model(model, act)
    fm, tm = CM.compliant_joint_loads_in_place!(F, T, cmodel, 0, pos, quat, vel, ω, mount, restv)
    dx = CM.compliant_multibody_dynamics(model, x, 0.0; joint_rest_quaternions=rest, joint_actuators=act)
    dpos = zeros(3, 3); dquat = zeros(4, 3); dvel = zeros(3, 3); dω = zeros(3, 3)
    CM.compliant_body_derivatives!(dpos, dquat, dvel, dω, cmodel, 0, F, T, quat, vel, ω)
    for i in 1:3
        b = 13 * (i - 1)
        @test dpos[:, i] ≈ dx[(b + 1):(b + 3)] rtol = 1e-14
        @test dquat[:, i] ≈ dx[(b + 4):(b + 7)] rtol = 1e-14
        @test dvel[:, i] ≈ dx[(b + 8):(b + 10)] rtol = 1e-13
        @test dω[:, i] ≈ dx[(b + 11):(b + 13)] rtol = 1e-13
    end
    # the reaction on the (fixed) mount is minus the force the first joint exerts on its child
    @test fm ≈ loads[1].translation_force_parent_world rtol = 1e-13
    @test (@allocated CM.compliant_joint_loads_in_place!(F, T, cmodel, 0, pos, quat, vel, ω, mount, restv)) == 0 skip=!_ALLOC_CHECKS
end

# ---------------------------------------------------------------------------
# (f) guards
# ---------------------------------------------------------------------------

@testset "guards" begin
    build = chain_build()
    sc, bus = rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b)])
    good = engine_args(sc; mission_time=1.0, solver_mode=:tsit5)
    @test SE._validate_attachments!(good, :tsit5) === nothing
    for mode in (:tsit5, :auto_stiff, :rodas5p, :dp8)
        @test SE._validate_attachments!(good, mode) === nothing
    end
    # constructor guards
    free = SM.build_compliant_topology([SM.CompliantTopologyNode(:a; mass_kg=1.0, inertia_body_kg_m2=0.1, position=Z3),
            SM.CompliantTopologyNode(:b; mass_kg=1.0, inertia_body_kg_m2=0.1, position=SVector(1.0, 0, 0))],
        [SM.CompliantTopologyEdge(:e, 1, 2; parent_point_body=Z3, child_point_body=Z3, k_translation_n_m=1.0)])
    @test_throws ArgumentError SM.CompliantAttachment(; model=free, link=bus)                       # no joint to the mount
    @test_throws ArgumentError SM.CompliantAttachment(; model=build, link=bus, initial_state=zeros(5))
    @test_throws ArgumentError SM.CompliantAttachment(; model=build, link=bus, rest_schedule=:not_callable)
    @test_throws ArgumentError SM.CompliantAttachment(; model=build, link=bus, mount_quaternion=zeros(4))
    @test_throws ArgumentError SM.CompliantAttachment(; model=build, link=bus, joint_actuators=[CM.CompliantJointActuator(:x, 99)])
    # link not in the spacecraft
    other = mklink(m=1.0, dims=(1.0, 1.0, 1.0))
    @test_throws ArgumentError SM.SpacecraftModel(; links=[bus], root=bus, attachments=[SM.CompliantAttachment(; model=build, link=other)])
    # ... also when the list is edited after construction
    bad_sc, _ = rigid_sc(SM.CompliantAttachment[])
    push!(bad_sc.attachments, SM.CompliantAttachment(; model=build, link=other))
    @test_throws ArgumentError SE._validate_attachments!(engine_args(bad_sc; mission_time=1.0), :dp8)
    # orientation_sim = false
    @test_throws ArgumentError SpaceAGORA.run_simulation(engine_args(sc; mission_time=1.0, orientation=false))
    # unsupported solver modes
    for mode in (:split_imex, :multirate, :symplectic, :gravity_backbone_split)
        @test_throws ArgumentError SE._validate_attachments!(good, mode)
    end
    @test_throws ArgumentError SpaceAGORA.run_simulation(engine_args(sc; mission_time=1.0, solver_mode=:split_imex))
    # forced flat route
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "flat") do
        @test_throws ArgumentError SpaceAGORA.run_simulation(engine_args(sc; mission_time=1.0))
    end
    # robot arm on the same spacecraft
    arm_model = SM.default_cloth_arm_model(link_lengths_m=(0.9, 0.8, 0.6), link_radii_m=(0.06, 0.05, 0.04), link_masses_kg=(6.0, 4.0, 2.0), mount_offset_body=(0.5, 0.0, 0.6))
    abase = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    target = SM.cloth_fk(arm_model, abase, [-0.12, -0.08, 0.06]).end_effector_position
    plan = SM.plan_robot_arm_motion(arm_model, abase, [0.08, 0.95, -0.85], target; config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=2.0))
    arm = SM.RobotArmControlEffector(plan=plan, spacecraft_idx=1, controller=SM.init_robot_arm_joint_mpc(plan; dt_s=0.1, horizon=6), control_dt_s=0.1)
    @test_throws ArgumentError SE._validate_attachments!(engine_args(sc; mission_time=1.0, control=(arm,), control_rates=[0.1]), :tsit5)
    # mixed runs: every spacecraft must carry the same attachment body count
    plain, _ = rigid_sc(SM.CompliantAttachment[])
    @test_throws ArgumentError SE._validate_attachments!(engine_args([sc, plain]; mission_time=1.0), :tsit5)
    two, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=chain_build(n=2), link=b)])
    @test_throws ArgumentError SE._validate_attachments!(engine_args([sc, two]; mission_time=1.0), :tsit5)
    same, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b)], ic=ic_at((0.0, 7.0e6, 0.0), (-7.5e3, 0.0, 0.0)))
    @test SE._validate_attachments!(engine_args([sc, same]; mission_time=1.0), :tsit5) === nothing
    # the flat route reroutes in auto mode
    p = SM.ODEParams(n_sats=1, args=good)
    SE._initialize_attachment_runtimes!(p)
    @test p.shared_buffers.attachments_present[]
    flat = (mode=:flat_constellation_effector_queue, allotment=4, scheduler=:dynamic, dominant_axis=:flat_effector, policy_applied=true,
        effector_decision=(use_threads=true, allotment=4, mode=:threaded, policy_applied=true))
    @test SE._articulated_plan_guard(flat, p).mode == :satellite_batch
end

# ---------------------------------------------------------------------------
# state layout and a short two-satellite constellation run
# ---------------------------------------------------------------------------

@testset "state layout, default initial state, constellation" begin
    build = chain_build(n=3)
    # explicit initial_state default: the build's state in the mount frame (identity base pose here)
    bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
    att = SM.CompliantAttachment(; model=build, link=bus)
    @test att.initial_state == build.initial_state
    @test SM.attachment_body_count(att) == 3
    # a bare model starts at its rest geometry with zero velocities
    rest = CM.compliant_rest_state(build.model)
    bare = SM.CompliantAttachment(; model=build.model, link=bus)
    @test bare.initial_state == rest
    @test rest[1:3] ≈ [0.5, 0.0, 0.0] atol = 1e-14
    @test rest[13 * 2 + 1:13 * 2 + 3] ≈ [2.5, 0.0, 0.0] atol = 1e-14
    # a build with a base pose is expressed relative to that base
    shifted = SM.build_compliant_topology(
        [SM.CompliantTopologyNode(:a; mass_kg=1.0, inertia_body_kg_m2=0.1, position=SVector(3.0, 1.0, 0.0), quaternion=rotq([0, 0, 1], π / 2))],
        [SM.CompliantTopologyEdge(:e, 0, 1; parent_point_body=SVector(1.0, 0.0, 0.0), child_point_body=Z3, k_translation_n_m=1.0)];
        base_position=SVector(2.0, 1.0, 0.0), base_quaternion=rotq([0, 0, 1], π / 2))
    x = CM.compliant_state_in_mount_frame(shifted.model, shifted.initial_state)
    @test x[1:3] ≈ [0.0, -1.0, 0.0] atol = 1e-14
    # state vector shape and initial state: absolute inertial pose and velocity from the mount motion
    ic = ic_at((7.0e6, 0.0, 0.0), (10.0, 20.0, 30.0); q=rotq([0, 0, 1], π / 2), ω=(0.0, 0.0, 0.1))
    sc, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=build.model, link=b, mount_point=(1.0, 0.0, 0.0))]; ic=ic)
    args = engine_args(sc; mission_time=1.0)
    u = SE.build_initial_conditions(args)
    v = u.sc[1]
    @test size(v.att_r) == (3, 3) && size(v.att_q) == (4, 3) && size(v.att_v) == (3, 3) && size(v.att_ω) == (3, 3)
    # bus rotated +90 deg about z: the mount at link-frame (1,0,0) sits at +y; body 1 at mount-frame (0.5,0,0) -> (0, 1.5, 0)
    # att_r and att_v are RELATIVE to the bus position and velocity (inertial axes)
    @test v.att_r[:, 1] ≈ [0.0, 1.5, 0.0] atol = 1e-12
    # velocity relative to the bus: omega x r (r from the bus origin)
    @test v.att_v[:, 1] ≈ cross([0, 0, 0.1], [0.0, 1.5, 0.0]) atol = 1e-12
    @test v.att_ω[:, 1] ≈ [0.0, 0.0, 0.1] atol = 1e-14
    # a two-spacecraft run: same dynamics per satellite, independent states
    ic2 = ic_at((0.0, 7.0e6, 0.0), (-7.5e3, 0.0, 0.0))
    sc1, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=chain_build(n=3), link=b, mount_point=(0.5, 0.0, 0.0))]; ic=ic_at((7.0e6, 0.0, 0.0), (0.0, 7.5e3, 0.0)))
    sc2, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=chain_build(n=3, deflect=false), link=b, mount_point=(0.5, 0.0, 0.0))]; ic=ic2)
    two = engine_args([sc1, sc2]; mission_time=5.0, data_rate=1.0, effectors=(SM.InverseSquaredGravityModel(),), solver_mode=:tsit5, tol=tight_tolerances(dt_max=0.05, rel=1e-10))
    for mode in ("auto", "serial", "satellite")
        r = withenv("SPACEAGORA_RHS_EXECUTION_MODE" => mode) do
            SpaceAGORA.run_simulation(two; return_results=true).table
        end
        @test "sc2_attachment_pose_21" in names(r)
        @test r.sc2_attachment_pose_2[end] != r.sc1_attachment_pose_2[end]
    end
end

# ---------------------------------------------------------------------------
# (g) allocations
# ---------------------------------------------------------------------------

function _setup_alloc(sc, eff)
    args = engine_args(sc; mission_time=10.0, effectors=(eff,), solver_mode=:tsit5)
    p = SM.ODEParams(n_sats=1, args=args)
    SE._initialize_runtime_env_config!(p)
    SE._initialize_articulated_runtimes!(p)
    SE._initialize_attachment_runtimes!(p)
    u = SE.build_initial_conditions(args)
    return p, u, zero(u)
end

_rigid_alloc(p, u, du, forces, torques, eff) = @allocated SE._apply_attachments_rigid!(du.sc[1], u.sc[1], p, 1, 0.0, forces, torques, (eff,))
_art_alloc(art, p, u, du, forces, torques, eff) = @allocated SE._assign_articulated_rhs!(du.sc[1], u.sc[1], art, p, 1, 0.0, forces, torques, 0.0, (eff,))
_dispatch_alloc(p, u, du) = @allocated SE._spacecraft_dynamics_dispatch!(du, u, p, 0.0)

# time-varying rest schedule that allocates nothing
struct SwingRest end
function (::SwingRest)(out::Vector{SVector{4, Float64}}, t::Float64)
    s = 0.1 * sin(2t)
    @inbounds out[2] = SVector(0.0, 0.0, sin(s / 2), cos(s / 2))
    return nothing
end

@testset "attachment RHS is allocation-free" begin
    build = chain_build(n=3, kx=200.0, cx=2.0, kr=20.0, cr=0.5)
    for (label, schedule) in (("no schedule", nothing), ("allocation-free schedule", SwingRest()))
        for eff in (SM.InverseSquaredGravityModel(), SM.InverseSquaredJ2GravityModel())
            sc, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b, mount_point=(0.5, 0.2, 0.1), rest_schedule=schedule)];
                ic=ic_at((7.0e6, 0.0, 0.0), (0.0, 7.5e3, 0.0); ω=(0.0, 0.0, 1e-3)))
            p, u, du = _setup_alloc(sc, eff)
            forces = MVector{3, Float64}(0.0, 0.0, 0.0); torques = MVector{3, Float64}(0.0, 0.0, 0.0)
            SE._apply_attachments_rigid!(du.sc[1], u.sc[1], p, 1, 0.0, forces, torques, (eff,))     # warm-up
            alloc = _rigid_alloc(p, u, du, forces, torques, eff)
            @info "rigid attachment RHS allocation" label effector = nameof(typeof(eff)) alloc
            @test alloc == 0 skip=!_ALLOC_CHECKS
            @test any(!iszero, du.sc[1].att_v)                                  # the attachment actually has a derivative
        end
        # whole RHS dispatch, point-mass gravity: the attachment adds nothing to what the same rigid spacecraft allocates without it
        eff = SM.InverseSquaredGravityModel()
        ic = ic_at((7.0e6, 0.0, 0.0), (0.0, 7.5e3, 0.0); ω=(0.0, 0.0, 1e-3))
        sc, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b, mount_point=(0.5, 0.2, 0.1), rest_schedule=schedule)]; ic=ic)
        base_sc, _ = rigid_sc(SM.CompliantAttachment[]; ic=ic)
        allocs = map((sc, base_sc)) do s
            p, u, du = _setup_alloc(s, eff)
            SE._spacecraft_dynamics_dispatch!(du, u, p, 0.0)
            SE._spacecraft_dynamics_dispatch!(du, u, p, 0.0)
            _dispatch_alloc(p, u, du)
        end
        @info "full RHS dispatch allocation, with attachment vs without" label with_attachment = allocs[1] without = allocs[2]
        @test allocs[1] == allocs[2]

        # articulated spacecraft, attachment on the moving link
        for eff in (SM.InverseSquaredGravityModel(), SM.InverseSquaredJ2GravityModel())
            sca, _, _ = hinge_sc((bus, panel) -> [SM.CompliantAttachment(; model=build, link=panel, mount_point=(0.0, 0.5, 0.0), rest_schedule=schedule)])
            p, u, du = _setup_alloc(sca, eff)
            art = p.shared_buffers.articulated_runtimes[1]
            forces = MVector{3, Float64}(0.1, 0.2, 0.3); torques = MVector{3, Float64}(0.01, 0.0, 0.02)
            SE._assign_articulated_rhs!(du.sc[1], u.sc[1], art, p, 1, 0.0, forces, torques, 0.0, (eff,))
            alloc = _art_alloc(art, p, u, du, forces, torques, eff)
            @info "articulated attachment RHS allocation" label effector = nameof(typeof(eff)) alloc
            @test alloc == 0 skip=!_ALLOC_CHECKS
        end
    end
    # runs with no attachments pay only the flag: the flag is false and the helper returns immediately
    plain, _ = rigid_sc(SM.CompliantAttachment[])
    p, u, du = _setup_alloc(plain, SM.InverseSquaredGravityModel())
    @test !p.shared_buffers.attachments_present[]
    forces = MVector{3, Float64}(0.0, 0.0, 0.0); torques = MVector{3, Float64}(0.0, 0.0, 0.0)
    SE._apply_attachments_rigid!(du.sc[1], u.sc[1], p, 1, 0.0, forces, torques, (SM.InverseSquaredGravityModel(),))
    @test _rigid_alloc(p, u, du, forces, torques, SM.InverseSquaredGravityModel()) == 0 skip=!_ALLOC_CHECKS
end

# ---------------------------------------------------------------------------
# checkpoint / copies
# ---------------------------------------------------------------------------

@testset "attachments survive deepcopy and serialization" begin
    using Serialization
    build = chain_build()
    sc, _ = rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b, mount_point=(0.5, 0.0, 0.0), rest_schedule=SwingRest())])
    for c in (deepcopy(sc), (io = IOBuffer(); serialize(io, sc); seekstart(io); deserialize(io)))
        @test c.attachments[1].link === c.root
        @test c.attachments[1].initial_state == sc.attachments[1].initial_state
        @test SE._validate_attachments!(engine_args(c; mission_time=1.0), :dp8) === nothing
    end
end

# Indexed control effector (bound to its spacecraft when per-spacecraft GNC is flattened): a constant inertial force.
struct PushControl <: SM.AbstractControlEffectorModel
    sat_idx::Int
    force::SVector{3, Float64}
end
SM.calcControlEffect!(::PushControl, u, p, t::Float64, i::Int64) = nothing
SM.calcControlForceTorque(m::PushControl, u::AbstractVector, p::SM.ODEParams, i::Int64, t::Float64) = i == m.sat_idx ? (m.force, SVector(0.0, 0.0, 0.0)) : (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
SM.calcControlMassFlowRate(::PushControl, u::AbstractVector, p::SM.ODEParams, i::Int64, t::Float64)::Float64 = 0.0
SpaceAGORA.bind_spacecraft(m::PushControl, sat_idx::Int) = PushControl(sat_idx, m.force)

@testset "constructors and GNC flattening keep attachments" begin
    build = chain_build()
    sc, bus = rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b)])
    @test length(SE._without_gnc(sc).attachments) == 1 && SE._without_gnc(sc).attachments[1] === sc.attachments[1]
    # the 11- and 14-argument positional constructors still work and declare no attachments
    pos11 = SM.SpacecraftModel(sc.joints, sc.links, sc.root, sc.instant_actuation, sc.dry_mass, sc.prop_mass, sc.inertia_tensor,
        sc.n_reaction_wheels, sc.n_thrusters, sc.initial_condition, sc.id)
    pos14 = SM.SpacecraftModel(sc.joints, sc.links, sc.root, sc.instant_actuation, sc.dry_mass, sc.prop_mass, sc.inertia_tensor,
        sc.n_reaction_wheels, sc.n_thrusters, sc.initial_condition, sc.id, sc.guidance, sc.navigation, sc.control)
    @test isempty(pos11.attachments) && isempty(pos14.attachments)
    pos15 = SM.SpacecraftModel(sc.joints, sc.links, sc.root, sc.instant_actuation, sc.dry_mass, sc.prop_mass, sc.inertia_tensor,
        sc.n_reaction_wheels, sc.n_thrusters, sc.initial_condition, sc.id, sc.guidance, sc.navigation, sc.control, sc.attachments)
    @test pos15.attachments == sc.attachments
    # a spacecraft that declares its own GNC and carries an attachment runs through the flattening path
    gnc_sc = SM.SpacecraftModel(; links=[bus], root=bus, inertia_tensor=bus.inertia, initial_condition=sc.initial_condition, attachments=sc.attachments,
        control=SM.ControlModel(control_effectors=(PushControl(0, SVector(0.0, 5.0, 0.0)),), control_rates=[1.0]))
    res = SpaceAGORA.run_simulation(engine_args(gnc_sc; mission_time=1.0, data_rate=0.5, solver_mode=:tsit5); return_results=true)
    @test "sc1_attachment_pose_1" in names(res.table)
    @test res.table.sc1_vel_2[end] > res.table.sc1_vel_2[1] + 0.01        # the bound effector pushed the bus (5 N on 20 kg for 1 s)
end

@testset "checkpoint resume continues a run with an attachment" begin
    build = chain_build(n=3, kx=200.0, cx=1.0, kr=20.0, cr=0.2)
    sc() = first(rigid_sc(b -> [SM.CompliantAttachment(; model=build, link=b, mount_point=(0.5, 0.0, 0.0), rest_schedule=SwingRest())];
        ic=ic_at((7.0e6, 0.0, 0.0), (30.0, -20.0, 10.0); ω=(0.01, 0.02, -0.01))))
    function with_settings(a; overrides...)
        st = a.simulation_settings
        names = fieldnames(typeof(st))
        vals = NamedTuple{names}(map(n -> getfield(st, n), names))
        return SM.SimConfig._with_configuration(a; simulation_settings=SM.SimulationSettings(; merge(vals, overrides)...))
    end
    dir = mktempdir()
    cfg(mission; kw...) = with_settings(engine_args(sc(); mission_time=mission, data_rate=1.0, results_directory=dir, tol=tight_tolerances(dt_max=0.05)); kw...)
    SpaceAGORA.run_simulation(cfg(5.0; checkpoint_enabled=true, checkpoint_interval_s=2.0); return_results=true)
    resumed = SpaceAGORA.run_simulation(cfg(10.0; checkpoint_enabled=true, checkpoint_interval_s=2.0, resume_from_checkpoint=true); return_results=true)
    straight = SpaceAGORA.run_simulation(cfg(10.0); return_results=true)
    cols = vcat(["sc1_attachment_pose_$k" for k in 1:21], ["sc1_pos_$k" for k in 1:3], ["sc1_q_$k" for k in 1:4])
    diff = maximum(abs(resumed.table[end, c] - straight.table[end, c]) / max(1.0, abs(straight.table[end, c])) for c in cols)
    @info "attachment resume vs straight, max relative final-state difference" diff
    @test resumed.table.time[end] ≈ 10.0
    @test resumed.table.time[1] >= 5.0 - 1e-9
    @test diff < 1e-8
end

end # module
