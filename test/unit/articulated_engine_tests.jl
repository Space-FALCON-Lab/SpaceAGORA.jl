module ArticulatedEngineTests

# Engine integration of the articulated-body library: run_simulation with spacecraft whose
# Joints are not :fixed. Library-level tests live in articulated_body_tests.jl.

using Test
using LinearAlgebra
using StaticArrays
using DataFrames
using Serialization
using SpaceAGORA
using SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const AB = SM.ArticulatedBody
const SE = SpaceAGORA.SimulationEngine

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

rotq(axis, ang) = (a = normalize(collect(Float64, axis)); SVector{4, Float64}(a[1] * sin(ang / 2), a[2] * sin(ang / 2), a[3] * sin(ang / 2), cos(ang / 2)))
Rmat(q) = AB._rotmat(SVector{4, Float64}(q))

function mklink(; root=false, m, dims, r=(0.0, 0.0, 0.0), q=(0.0, 0.0, 0.0, 1.0), kwargs...)
    return SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(r...), q=MVector{4, Float64}(q...); kwargs...)
end

# Joint at point P (bus frame) between a parent and a child link at their configured geometry.
function mkjoint(par::SM.Link, child::SM.Link, P; kwargs...)
    rp = par.root ? zeros(3) : collect(par.r)
    Rp = par.root ? Matrix(1.0I, 3, 3) : Matrix(Rmat(par.q))
    Rc = Matrix(Rmat(child.q))
    p1 = SVector{3, Float64}(Rp' * (collect(P) - rp))
    p2 = SVector{3, Float64}(Rc' * (collect(P) - collect(child.r)))
    return SM.Joint(par, p1, child, p2; kwargs...)
end

const T0 = SM.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)

function tight_tolerances(; dt_max=0.05)
    return SM.IntegrationTolerances(
        reltol_orbit=1e-12, abstol_orbit=1e-14, reltol_atmosphere=1e-12, abstol_atmosphere=1e-14,
        reltol_quaternion=1e-12, abstol_quaternion=1e-14, reltol_mass=1e-12, abstol_mass=1e-14,
        reltol_angular_rate=1e-12, abstol_angular_rate=1e-14, dt_max_orbit=dt_max, dt_max_atmosphere=dt_max,
    )
end

function engine_args(sc; mission_time, data_rate=0.05, effectors=(), control=(), control_rates=Float64[],
        orientation=true, solver_mode=:tsit5, tol=tight_tolerances(), results_directory=mktempdir(), ephem=SM.SimpleEphemeridesModel())
    base = make_example_config(
        planet=SM.make_no_gram_planet(:earth), spacecraft=sc, mission_time=mission_time, initial_time=T0,
        dynamic_effectors=effectors, density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
        orientation_sim=orientation, keplerian=true, verbose=false, results=false, results_directory=results_directory,
        solver_config=SM.SolverConfig(solver_mode=solver_mode),
    )
    return SM.SimConfig._with_configuration(base;
        integration_tolerances=tol,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=orientation, num_steps_to_save=1000, data_rate=data_rate),
        control_model=SM.ControlModel(control_effectors=control, control_rates=control_rates),
    )
end

const IC0 = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(30.0, -20.0, 10.0);
    q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, 0.0))

# bus + one hinged panel (the stage-A free-flight case)
function hinge_spacecraft(; k, c, θ0, ic=IC0, prop=0.0, joint_type=:hinge)
    bus = mklink(root = true, m = 10.0, dims = (1.0, 1.0, 1.0))
    panel = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, 1.1, 0.0))
    joint = mkjoint(bus, panel, (0.0, 0.5, 0.0); joint_type = joint_type, axis = [0, 0, 1], stiffness = k, damping = c, initial_q = θ0)
    return SM.SpacecraftModel(; joints = [joint], links = [bus, panel], root = bus, initial_condition = ic, prop_mass = prop)
end

function hinge_ieff()
    Ib = 10.0 / 12 * 2.0
    Ip = 2.0 / 12 * (0.05^2 + 1.0^2)
    a, b = 0.5, 0.6
    μ = 10.0 * 2.0 / 12.0
    M11 = Ib + Ip + μ * (a + b)^2; M12 = Ip + μ * b * (a + b); M22 = Ip + μ * b^2
    return M22 - M12^2 / M11
end

function crossings(ts, ys)
    out = Float64[]
    for i in 2:length(ts)
        if ys[i - 1] * ys[i] < 0
            push!(out, ts[i - 1] - ys[i - 1] * (ts[i] - ts[i - 1]) / (ys[i] - ys[i - 1]))
        end
    end
    return out
end

# The stage-A chain: hinge + slide + ball with two fixed links merged in.
function chain_spacecraft(; ic=IC0, prop=0.5, initial_q=[0.4, 0.15], initial_qd=[0.5, 0.3], ball_q=rotq([1, 2, 1], 0.5), ball_qd=[0.2, -0.1, 0.15],
        k=(6.0, 15.0, 3.0), c=(0.0, 0.0, 0.0))
    bus = mklink(root = true, m = 20.0, dims = (1.0, 0.8, 0.6))
    A = mklink(m = 3.0, dims = (0.1, 1.2, 0.6), r = (0.0, 1.0, 0.2), q = rotq([1, 0, 0], 0.3))
    B = mklink(m = 2.0, dims = (0.2, 0.9, 0.4), r = (0.3, 2.2, 0.3), q = rotq([0, 1, 0], -0.5))
    C = mklink(m = 1.5, dims = (0.3, 0.3, 0.6), r = (0.4, 3.0, 0.1), q = rotq([1, 1, 0], 0.8))
    D = mklink(m = 2.0, dims = (0.4, 0.4, 0.4), r = (-0.8, 0.0, 0.1), q = rotq([0, 0, 1], 0.2))
    E = mklink(m = 0.8, dims = (0.1, 0.5, 0.3), r = (0.5, 1.2, -0.3), q = rotq([1, 1, 1], -0.6))
    joints = [
        mkjoint(bus, A, (0.0, 0.6, 0.1); joint_type = :hinge, axis = [0.3, 0.2, 0.9], stiffness = k[1], damping = c[1], initial_q = initial_q[1], initial_qd = initial_qd[1]),
        mkjoint(A, B, (0.1, 1.7, 0.25); joint_type = :slide, axis = [1.0, 0.5, 0.0], stiffness = k[2], damping = c[2], initial_q = initial_q[2], initial_qd = initial_qd[2]),
        mkjoint(B, C, (0.35, 2.7, 0.2); joint_type = :ball, stiffness = k[3], damping = c[3], rest = rotq([0, 0, 1], 0.2), initial_q = ball_q, initial_qd = ball_qd),
        mkjoint(bus, D, (-0.5, 0.0, 0.05)),
        mkjoint(A, E, (0.3, 1.0, -0.1)),
    ]
    return SM.SpacecraftModel(; joints = joints, links = [bus, A, B, C, D, E], root = bus, initial_condition = ic, prop_mass = prop)
end

# ---------------------------------------------------------------------------
# (a) all-fixed joints change nothing
# ---------------------------------------------------------------------------

@testset "all-fixed joints leave the run unchanged" begin
    function build(with_joints)
        bus = mklink(root = true, m = 20.0, dims = (1.0, 0.8, 0.6))
        L = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, -1.1, 0.0))
        R = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, 1.1, 0.0))
        joints = with_joints ? [mkjoint(bus, L, (0.0, -0.5, 0.0)), mkjoint(bus, R, (0.0, 0.5, 0.0))] : SM.Joint[]
        ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 7.5e3, 0.0); q = SVector(0.0, 0.0, 0.0, 1.0), ang_vel = SVector(0.0, 0.0, 1e-3))
        return SM.SpacecraftModel(; joints = joints, links = [bus, L, R], root = bus, initial_condition = ic, inertia_tensor = bus.inertia, prop_mass = 1.0)
    end
    run_it(sc) = SM.SimConfig._with_configuration(engine_args(sc; mission_time = 60.0, data_rate = 5.0, effectors = (SM.InverseSquaredJ2GravityModel(),),
        tol = SM.IntegrationTolerances()))
    with = SpaceAGORA.run_simulation(run_it(build(true)); return_results = true)
    without = SpaceAGORA.run_simulation(run_it(build(false)); return_results = true)
    @test names(with.table) == names(without.table)
    @test all(isequal(with.table[!, n], without.table[!, n]) for n in names(with.table))
    @test !any(startswith("sc1_joint"), names(with.table))
    args = run_it(build(true))
    p = SM.ODEParams(n_sats = 1, args = args)
    SE._initialize_articulated_runtimes!(p)
    @test all(isnothing, p.shared_buffers.articulated_runtimes)
    @test !p.shared_buffers.articulated_present[]
    @test !SE._spacecraft_articulated(args.dynamics_model.spacecraft[1])
end

# ---------------------------------------------------------------------------
# (b) free-flight hinge through run_simulation
# ---------------------------------------------------------------------------

@testset "engine free-flight hinge: frequency and decay" begin
    k = 4.0
    Ieff = hinge_ieff()
    Tper = 2π / sqrt(k / Ieff)
    tend = 12Tper
    res = SpaceAGORA.run_simulation(engine_args(hinge_spacecraft(k = k, c = 0.0, θ0 = 0.01); mission_time = tend, data_rate = 0.02, solver_mode = :dp8); return_results = true)
    df = res.table
    zc = crossings(df.time, df.sc1_joint_q_1)
    n = length(zc)
    @test n >= 20
    Tmeas = 2 * (zc[end] - zc[1]) / (n - 1)
    ferr = abs(Tmeas - Tper) / Tper
    @info "engine hinge frequency" Tpred = Tper Tmeas ferr
    @test ferr < 0.005

    c = 0.1 * 2 * sqrt(k * Ieff)
    λ = c / (2Ieff)
    resd = SpaceAGORA.run_simulation(engine_args(hinge_spacecraft(k = k, c = c, θ0 = 0.01); mission_time = tend, data_rate = 0.02, solver_mode = :dp8); return_results = true)
    dfd = resd.table
    te = crossings(dfd.time, dfd.sc1_joint_qd_1)
    amp = Float64[]
    for t in te
        i = searchsortedlast(dfd.time, t)
        w = (t - dfd.time[i]) / (dfd.time[i + 1] - dfd.time[i])
        push!(amp, abs((1 - w) * dfd.sc1_joint_q_1[i] + w * dfd.sc1_joint_q_1[i + 1]))
    end
    slope = ([ones(length(te)) te] \ log.(amp))[2]
    derr = abs(-slope - λ) / λ
    @info "engine hinge decay" λpred = λ λmeas = -slope derr
    @test derr < 0.02
    # Saved outputs: joint state, link poses and system COM columns exist with the documented sizes.
    @test all(n -> n in names(dfd), ["sc1_joint_q_1", "sc1_joint_qd_1", "sc1_system_com_3"])
    @test count(startswith("sc1_articulated_link_pose_"), names(dfd)) == 14
    # Free flight: the system COM does not accelerate, so it moves at the initial velocity.
    v0 = [30.0, -20.0, 10.0]
    com_rate = [(dfd[end, "sc1_system_com_$i"] - dfd[1, "sc1_system_com_$i"]) / (dfd.time[end] - dfd.time[1]) for i in 1:3]
    @test norm(com_rate - v0) < 1e-3
    @test dfd.sc1_mass[end] == dfd.sc1_mass[1]
end

# ---------------------------------------------------------------------------
# (c) conservation with a hinge + slide + ball chain
# ---------------------------------------------------------------------------

function sys_invariants(tree, art_moving_mass, sc_state)
    base = AB.ArticulatedBaseState(zeros(3), sc_state.vel, sc_state.q, sc_state.ω)   # translation invariant: origin at 0
    kin = AB.articulated_kinematics(tree, base, sc_state.joint_q, sc_state.joint_qd)
    masses = copy(tree.mass); masses[1] = sc_state.mass - art_moving_mass
    mt = sum(masses)
    X = sum(masses[b] * kin.pos[b] for b in 1:tree.nb) / mt
    V = sum(masses[b] * kin.vel[b] for b in 1:tree.nb) / mt
    L = zero(SVector{3, Float64}); KE = 0.0
    for b in 1:tree.nb
        R = Rmat(kin.quat[b]); ωb = kin.ω[b]
        L += masses[b] * cross(kin.pos[b] - X, kin.vel[b] - V) + R * (tree.inertia[b] * ωb)
        KE += masses[b] * dot(kin.vel[b], kin.vel[b]) / 2 + dot(ωb, tree.inertia[b] * ωb) / 2
    end
    return mt * V, L, KE + AB.articulated_potential_energy(tree, sc_state.joint_q)
end

@testset "engine conservation, hinge+slide+ball chain" begin
    sc = chain_spacecraft(ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(30.0, -20.0, 10.0); q = rotq([1, -1, 2], 0.9), ang_vel = SVector(0.05, 0.1, -0.04)))
    tree = AB.build_articulated_tree(sc)
    moving = AB.articulated_moving_mass(tree)
    args = engine_args(sc; mission_time = 100.0, data_rate = 1.0, solver_mode = :dp8, tol = tight_tolerances(dt_max = 0.1))
    res = SpaceAGORA.run_simulation(args; return_results = true, return_solution = true)
    sol = res.solution
    P0, L0, E0 = sys_invariants(tree, moving, sol.u[1].sc[1])
    dP = 0.0; dL = 0.0; dE = 0.0
    for u in sol.u
        P, L, E = sys_invariants(tree, moving, u.sc[1])
        dP = max(dP, norm(P - P0) / norm(P0)); dL = max(dL, norm(L - L0) / norm(L0)); dE = max(dE, abs(E - E0) / abs(E0))
    end
    @info "engine conservation drifts over 100 s" linear_momentum = dP angular_momentum = dL energy = dE steps = length(sol.t)
    @test dP < 1e-9
    @test dL < 1e-9
    @test dE < 1e-9
    @test abs(sol.u[end].sc[1].joint_q[1] - 0.4) > 1e-2      # the joints actually moved
    # Ball quaternion stays unit (projection callback).
    @test abs(norm(sol.u[end].sc[1].joint_q[3:6]) - 1) < 1e-12
    # Saved link poses: 6 links x (position, quaternion); link 1 is the root.
    row = res.table[end, :]
    pose = [row["sc1_articulated_link_pose_$k"] for k in 1:42]
    @test all(l -> abs(norm(pose[(7l - 3):(7l)]) - 1) < 1e-9, 1:6)
    # Link 1 is the root link; the state `pos` is the root composite COM (the fixed link D shifts it).
    @test norm(pose[1:3] - [row.sc1_pos_1, row.sc1_pos_2, row.sc1_pos_3]) ≈ norm(tree.link_com_in_body[1]) atol = 1e-7
    # The RHS never mutates the configured Link r/q.
    before = [(copy(l.r), copy(l.q)) for l in sc.links]
    SpaceAGORA.run_simulation(engine_args(sc; mission_time = 5.0, data_rate = 1.0); isolate_state = false)
    @test all(i -> sc.links[i].r == before[i][1] && sc.links[i].q == before[i][2], eachindex(sc.links))
end

# ---------------------------------------------------------------------------
# (d) circular orbit, stiff damped hinge vs. the same spacecraft all-fixed
# ---------------------------------------------------------------------------

@testset "orbit: stiff hinge vs. all-fixed" begin
    planet = SM.make_no_gram_planet(:earth)
    r0 = 7.0e6
    vc = sqrt(planet.μ / r0)
    n = vc / r0
    ic = SM.CartesianInitialCondition(SVector(r0, 0.0, 0.0), SVector(0.0, vc, 0.0); q = SVector(0.0, 0.0, 0.0, 1.0), ang_vel = SVector(0.0, 0.0, n))
    function three_body(; moving)
        bus = mklink(root = true, m = 20.0, dims = (1.0, 0.8, 0.6))
        L = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, -1.1, 0.0))
        R = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, 1.1, 0.0))
        kw = moving ? (joint_type = :hinge, axis = [0, 0, 1], stiffness = 400.0, damping = 40.0) : NamedTuple()
        joints = [mkjoint(bus, L, (0.0, -0.5, 0.0); kw...), mkjoint(bus, R, (0.0, 0.5, 0.0); kw...)]
        return bus, SM.SpacecraftModel(; joints = joints, links = [bus, L, R], root = bus, initial_condition = ic)
    end
    _, fixed_sc = three_body(moving = false)
    composite = AB.build_articulated_tree(fixed_sc).inertia[1]       # all-fixed composite about the system COM (bus origin by symmetry)
    fixed_sc.inertia_tensor = composite
    _, art_sc = three_body(moving = true)
    T = 2π / n
    tend = T                                                          # one orbit
    run_it(sc) = SpaceAGORA.run_simulation(engine_args(sc; mission_time = tend, data_rate = 20.0, effectors = (SM.InverseSquaredGravityModel(),),
        solver_mode = :dp8, tol = tight_tolerances(dt_max = 5.0)); return_results = true)
    art = run_it(art_sc)
    rig = run_it(fixed_sc)
    @test nrow(art.table) == nrow(rig.table)
    dcom = maximum(norm([art.table[i, "sc1_system_com_$k"] - rig.table[i, "sc1_pos_$k"] for k in 1:3]) for i in 1:nrow(art.table))
    dq = maximum(norm([art.table[i, "sc1_q_$k"] - rig.table[i, "sc1_q_$k"] for k in 1:4]) for i in 1:nrow(art.table))
    θmax = maximum(abs, art.table.sc1_joint_q_1)
    @info "orbit stiff hinge vs all-fixed" system_com_max_diff_m = dcom bus_quaternion_max_diff = dq max_hinge_angle_rad = θmax
    # The aligned LVLH configuration is an equilibrium of both models; the residual is the tidal torque on the hinge
    # over the hinge stiffness plus the second-order (d/r)^2 difference of point-mass gravity vs. the rigid model.
    @test θmax < 1e-10
    @test dcom < 1e-4
    @test dq < 1e-6
end

@testset "orbit: misaligned attitude, isotropic bodies, vs. rigid with gravity gradient" begin
    # Every body is a cube (isotropic inertia), so the intrinsic gravity-gradient torque the articulated model
    # omits is exactly zero; the rigid comparator carries the composite inertia and the gravity-gradient flag.
    planet = SM.make_no_gram_planet(:earth)
    r0 = 7.0e6
    vc = sqrt(planet.μ / r0)
    n = vc / r0
    ic = SM.CartesianInitialCondition(SVector(r0, 0.0, 0.0), SVector(0.0, vc, 0.0); q = rotq([1, 2, 3], 0.3), ang_vel = SVector(1e-3, 2e-3, n))
    function cubes(; moving)
        bus = mklink(root = true, m = 20.0, dims = (1.0, 1.0, 1.0))
        L = mklink(m = 2.0, dims = (0.2, 0.2, 0.2), r = (0.0, -1.1, 0.0))
        R = mklink(m = 2.0, dims = (0.2, 0.2, 0.2), r = (0.0, 1.1, 0.0))
        kw = moving ? (joint_type = :hinge, axis = [0, 0, 1], stiffness = 400.0, damping = 40.0) : NamedTuple()
        joints = [mkjoint(bus, L, (0.0, -0.5, 0.0); kw...), mkjoint(bus, R, (0.0, 0.5, 0.0); kw...)]
        return SM.SpacecraftModel(; joints = joints, links = [bus, L, R], root = bus, initial_condition = ic)
    end
    fixed_sc = cubes(moving = false)
    fixed_sc.inertia_tensor = AB.build_articulated_tree(fixed_sc).inertia[1]
    art_sc = cubes(moving = true)
    T = 2π / n
    run_it(sc, eff) = SpaceAGORA.run_simulation(engine_args(sc; mission_time = T, data_rate = 20.0, effectors = (eff,),
        solver_mode = :dp8, tol = tight_tolerances(dt_max = 5.0)); return_results = true)
    art = run_it(art_sc, SM.InverseSquaredGravityModel())
    rig = run_it(fixed_sc, SM.InverseSquaredGravityModel(gravity_gradient = true))
    dcom = maximum(norm([art.table[i, "sc1_system_com_$k"] - rig.table[i, "sc1_pos_$k"] for k in 1:3]) for i in 1:nrow(art.table))
    dq = maximum(norm([art.table[i, "sc1_q_$k"] - rig.table[i, "sc1_q_$k"] for k in 1:4]) for i in 1:nrow(art.table))
    θmax = maximum(abs, art.table.sc1_joint_q_1)
    @info "misaligned cubes: stiff hinge vs rigid + gravity gradient" system_com_max_diff_m = dcom bus_quaternion_max_diff = dq max_hinge_angle_rad = θmax
    # The gravity-gradient torque twists the attitude by O(1) over the orbit in both models, so the agreement is a real test.
    @test maximum(abs, rig.table.sc1_q_1 .- rig.table.sc1_q_1[1]) > 1e-2
    @test θmax < 1e-7                    # tidal torque on the hinge / stiffness
    @test dcom < 1e-4                    # dominated by the 1e-12 relative tolerance on a 7e6 m position
    @test dq < 1e-7
end

# ---------------------------------------------------------------------------
# (e) mass flow
# ---------------------------------------------------------------------------

struct ConstantThrust <: SM.AbstractControlEffectorModel
    force::SVector{3, Float64}
    mdot::Float64
end
SM.calcControlEffect!(::ConstantThrust, u, p, t::Float64, i::Int64) = nothing
SM.calcControlForceTorque(m::ConstantThrust, u::AbstractVector, p::SM.ODEParams, i::Int64, t::Float64) = (m.force, SVector(0.0, 0.0, 0.0))
SM.calcControlMassFlowRate(m::ConstantThrust, u::AbstractVector, p::SM.ODEParams, i::Int64, t::Float64)::Float64 = m.mdot

@testset "mass flow lowers the root mass only" begin
    sc = hinge_spacecraft(k = 4.0, c = 0.0, θ0 = 0.0, prop = 5.0)
    thrust = ConstantThrust(SVector(0.5, 0.0, 0.0), -1e-2)
    res = SpaceAGORA.run_simulation(engine_args(sc; mission_time = 100.0, data_rate = 10.0, control = (thrust,), control_rates = [1.0]); return_results = true)
    m0 = 10.0 + 2.0 + 5.0
    @test res.table.sc1_mass[1] ≈ m0
    @test res.table.sc1_mass[end] ≈ m0 - 1.0 rtol = 1e-9
    tree = AB.build_articulated_tree(sc)
    @test AB.articulated_moving_mass(tree) == 2.0               # the panel mass never changes
    # Runtime root-mass argument: depleting 3 kg is the same as building the tree with 3 kg less propellant.
    ws = AB.ArticulatedWorkspace(tree)
    base = AB.ArticulatedBaseState(zeros(3), zeros(3), [0, 0, 0, 1.0], [0.1, 0.2, 0.3])
    g = r -> SVector(0.0, 0.0, -1.0)
    F = SVector(1.0, 2.0, 3.0); Tq = SVector(0.1, 0.0, 0.2)
    a1, α1, q1 = AB.articulated_dynamics!(ws, tree, base, [0.2], [0.1], F, Tq, g; root_mass = tree.mass[1] - 3.0)
    a1c, α1c, q1c = collect(a1), collect(α1), copy(q1)
    tree2 = AB.build_articulated_tree(sc; prop_mass = 2.0)
    a2, α2, q2 = AB.articulated_dynamics!(AB.ArticulatedWorkspace(tree2), tree2, base, [0.2], [0.1], F, Tq, g)
    @test a1c ≈ collect(a2) rtol = 1e-14
    @test α1c ≈ collect(α2) rtol = 1e-14
    @test q1c ≈ q2 rtol = 1e-14
end

# ---------------------------------------------------------------------------
# (f) guards
# ---------------------------------------------------------------------------

function _arm_plan(; duration_s=10.0)
    arm = SM.default_cloth_arm_model(link_lengths_m=(0.9, 0.8, 0.6), link_radii_m=(0.06, 0.05, 0.04), link_masses_kg=(6.0, 4.0, 2.0), mount_offset_body=(0.5, 0.0, 0.6))
    base = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    target = SM.cloth_fk(arm, base, [-0.12, -0.08, 0.06]).end_effector_position
    return SM.plan_robot_arm_motion(arm, base, [0.08, 0.95, -0.85], target; config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=duration_s))
end

@testset "guards" begin
    sc = hinge_spacecraft(k = 4.0, c = 0.0, θ0 = 0.01)
    good = engine_args(sc; mission_time = 1.0, data_rate = 0.5)
    @test SE._validate_articulated_spacecraft!(good, :tsit5) === nothing

    # orientation_sim=false
    bad = engine_args(sc; mission_time = 1.0, orientation = false)
    @test_throws ArgumentError SpaceAGORA.run_simulation(bad)
    # unsupported solver modes
    for mode in (:split_imex, :multirate, :symplectic, :gravity_backbone_split)
        @test_throws ArgumentError SE._validate_articulated_spacecraft!(good, mode)
    end
    @test_throws ArgumentError SpaceAGORA.run_simulation(engine_args(sc; mission_time = 1.0, solver_mode = :split_imex))
    for mode in (:tsit5, :auto_stiff, :rodas5p, :dp8)
        @test SE._validate_articulated_spacecraft!(good, mode) === nothing
    end
    # forced flat route
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "flat") do
        @test_throws ArgumentError SpaceAGORA.run_simulation(engine_args(sc; mission_time = 1.0))
    end
    # reaction wheels
    rw_bus = mklink(root = true, m = 10.0, dims = (1.0, 1.0, 1.0), rw = [0.0], J_rw = reshape([1.0, 0.0, 0.0], 3, 1))
    panel = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, 1.1, 0.0))
    rw_sc = SM.SpacecraftModel(; joints = [mkjoint(rw_bus, panel, (0.0, 0.5, 0.0); joint_type = :hinge, axis = [0, 0, 1])], links = [rw_bus, panel], root = rw_bus, initial_condition = IC0)
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(rw_sc; mission_time = 1.0), :tsit5)
    # robot arm on the same spacecraft
    plan = _arm_plan(duration_s = 2.0)
    arm = SM.RobotArmControlEffector(plan = plan, spacecraft_idx = 1, controller = SM.init_robot_arm_joint_mpc(plan; dt_s = 0.1, horizon = 6), control_dt_s = 0.1)
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(sc; mission_time = 1.0, control = (arm,), control_rates = [0.1]), :tsit5)
    # dry_mass must equal the sum of the link masses
    off = hinge_spacecraft(k = 4.0, c = 0.0, θ0 = 0.0)
    off.dry_mass = 99.0
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(off; mission_time = 1.0), :tsit5)
    # kinematic panel-angle control: allowed on a link of the root body, an error on a moving link
    bus = mklink(root = true, m = 10.0, dims = (1.0, 1.0, 1.0))
    fixed_panel = mklink(m = 1.0, dims = (0.05, 1.0, 0.5), r = (0.0, -1.1, 0.0))
    moving_panel = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, 1.1, 0.0))
    three = SM.SpacecraftModel(; joints = [mkjoint(bus, fixed_panel, (0.0, -0.5, 0.0)),
            mkjoint(bus, moving_panel, (0.0, 0.5, 0.0); joint_type = :hinge, axis = [0, 0, 1], stiffness = 1.0)],
        links = [bus, fixed_panel, moving_panel], root = bus, initial_condition = IC0)
    ctl(links) = (SM.SolarPanelAngleOfAttackControlModel(controlled_panel_links = links),)
    @test SE._validate_articulated_spacecraft!(engine_args(three; mission_time = 1.0, control = ctl((2,)), control_rates = [1.0]), :tsit5) === nothing
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(three; mission_time = 1.0, control = ctl((3,)), control_rates = [1.0]), :tsit5)
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(three; mission_time = 1.0, control = ctl((2, 3)), control_rates = [1.0]), :tsit5)
    # Tree validation errors surface at setup too (mismatched attachment point on a moving joint).
    bad_joint = SM.Joint(bus, SVector(0.0, 0.5, 0.0), moving_panel, SVector(0.0, -0.6 - 1e-3, 0.0); joint_type = :hinge, axis = [0, 0, 1])
    bad_sc = SM.SpacecraftModel(; joints = [bad_joint], links = [bus, moving_panel], root = bus, initial_condition = IC0)
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(bad_sc; mission_time = 1.0), :tsit5)
end

@testset "flat queue reroutes, state survives deepcopy and serialization" begin
    sc = hinge_spacecraft(k = 4.0, c = 0.0, θ0 = 0.01)
    args = engine_args(sc; mission_time = 1.0)
    p = SM.ODEParams(n_sats = 1, args = args)
    SE._initialize_articulated_runtimes!(p)
    @test p.shared_buffers.articulated_present[]
    flat = (mode = :flat_constellation_effector_queue, allotment = 4, scheduler = :dynamic, dominant_axis = :flat_effector, policy_applied = true,
        effector_decision = (use_threads = true, allotment = 4, mode = :threaded, policy_applied = true))
    @test SE._articulated_plan_guard(flat, p).mode == :satellite_batch
    @test SE._articulated_plan_guard(flat, SM.ODEParams(n_sats = 1, args = args)).mode == :flat_constellation_effector_queue
    # Joint configuration survives deepcopy (the run isolates state) and process-route serialization.
    for copy_sc in (deepcopy(sc), (io = IOBuffer(); serialize(io, sc); seekstart(io); deserialize(io)))
        j = copy_sc.joints[1]
        @test j.joint_type === :hinge && j.stiffness == 4.0 && j.initial_q == 0.01
        @test j.link1 === copy_sc.links[findfirst(l -> l.root, copy_sc.links)]
        t1 = AB.build_articulated_tree(copy_sc)
        t0 = AB.build_articulated_tree(sc)
        @test t1.mass == t0.mass && t1.inertia == t0.inertia && t1.q0 == t0.q0
    end
end

@testset "articulated constellation, mixed-run guard" begin
    sc1 = hinge_spacecraft(k = 4.0, c = 0.2, θ0 = 0.05)
    sc2 = hinge_spacecraft(k = 6.0, c = 0.1, θ0 = -0.03, ic = SM.CartesianInitialCondition(SVector(0.0, 7.0e6, 0.0), SVector(-7.5e3, 0.0, 0.0); q = SVector(0.0, 0.0, 0.0, 1.0), ang_vel = SVector(0.0, 0.0, 0.0)))
    one = engine_args(sc1; mission_time = 30.0, data_rate = 5.0, effectors = (SM.InverseSquaredJ2GravityModel(),), tol = SM.IntegrationTolerances(reltol_orbit = 1e-10, abstol_orbit = 1e-12, reltol_angular_rate = 1e-10, abstol_angular_rate = 1e-12, reltol_quaternion = 1e-10, abstol_quaternion = 1e-12, dt_max_orbit = 1.0))
    two = SM.SimConfig._with_configuration(one; dynamics_model = SM.DynamicsModel([sc1, sc2], (SM.InverseSquaredJ2GravityModel(),)))
    r1 = SpaceAGORA.run_simulation(one; return_results = true).table
    for mode in ("auto", "serial", "satellite")
        r2 = withenv("SPACEAGORA_RHS_EXECUTION_MODE" => mode) do
            SpaceAGORA.run_simulation(two; return_results = true).table
        end
        @test "sc2_joint_q_1" in names(r2)
        @test maximum(abs.(r2.sc1_joint_q_1 .- r1.sc1_joint_q_1)) < 1e-7      # same dynamics, steps coupled through the shared solver
        @test r2.sc2_joint_q_1[end] != r2.sc1_joint_q_1[end]
    end
    # An articulated spacecraft cannot share a run with a rigid one (equal-sized state blocks).
    rigid_bus = mklink(root = true, m = 12.0, dims = (1.0, 1.0, 1.0))
    rigid = SM.SpacecraftModel(; joints = SM.Joint[], links = [rigid_bus], root = rigid_bus, initial_condition = IC0, inertia_tensor = rigid_bus.inertia)
    mixed = SM.SimConfig._with_configuration(one; dynamics_model = SM.DynamicsModel([sc1, rigid], (SM.InverseSquaredJ2GravityModel(),)))
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(mixed, :tsit5)
end

# ---------------------------------------------------------------------------
# base torque is about the root composite COM, not the bus origin
# ---------------------------------------------------------------------------

function _base_rhs(sc, F, torque_origin)
    args = engine_args(sc; mission_time = 10.0)
    p = SM.ODEParams(n_sats = 1, args = args)
    SE._initialize_runtime_env_config!(p)
    SE._initialize_articulated_runtimes!(p)
    u = SE.build_initial_conditions(args)
    du = zero(u)
    art = p.shared_buffers.articulated_runtimes[1]
    SE._assign_articulated_rhs!(du.sc[1], u.sc[1], art, p, 1, 0.0, MVector{3, Float64}(F), MVector{3, Float64}(torque_origin), 0.0, ())
    sv = u.sc[1]
    base = AB.ArticulatedBaseState(collect(sv.pos), collect(sv.vel), collect(sv.q), collect(sv.ω))
    a, α, qdd = AB.articulated_dynamics!(art.ws, art.tree, base, sv.joint_q, sv.joint_qd, F, zeros(3), r -> zeros(3); root_mass = sv.mass - art.moving_mass)
    return (tree = art.tree, du = du.sc[1], a = collect(a), α = collect(α), qdd = copy(qdd))
end

@testset "articulated base torque is moved from the bus origin to the root composite COM" begin
    ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(30.0, -20.0, 10.0);
        q = rotq([1, 1, 0], 0.7), ang_vel = SVector(0.0, 0.0, 0.0))
    bus = mklink(root = true, m = 10.0, dims = (1.0, 1.0, 1.0))
    block = mklink(m = 5.0, dims = (0.4, 0.4, 0.4), r = (0.5, 0.0, 0.0))
    panel = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, 1.1, 0.0))
    fixed = mkjoint(bus, block, (0.5, 0.0, 0.0); joint_type = :fixed)
    hinge = mkjoint(bus, panel, (0.0, 0.5, 0.0); joint_type = :hinge, axis = [0, 0, 1], stiffness = 1.0, damping = 0.0, initial_q = 0.0)
    sc = SM.SpacecraftModel(; joints = [fixed, hinge], links = [bus, block, panel], root = bus, initial_condition = ic)
    F = SVector(1.0, -2.0, 3.0)
    # a force through the composite COM: its torque about the bus origin is c x F_body
    R = Rmat(ic.q)
    tree0 = AB.build_articulated_tree(sc)
    c = tree0.root_com_bus
    @test norm(c) > 0.1
    r = _base_rhs(sc, F, cross(c, R' * F))
    @test collect(r.du.ω) ≈ r.α rtol = 1e-12
    @test collect(r.du.vel) ≈ r.a rtol = 1e-12
    @test collect(r.du.joint_qd) ≈ r.qdd rtol = 1e-12
    # Without the conversion the same load would spin the base: the unconverted torque is not neutral.
    r0 = _base_rhs(sc, F, zeros(3))
    @test norm(collect(r0.du.ω) - r.α) > 1e-6

    # No :fixed child mass: root_com_bus = 0 and the torque passes through bit for bit.
    sc2 = hinge_spacecraft(k = 4.0, c = 0.0, θ0 = 0.1, ic = ic)
    Tq = SVector(0.01, 0.02, -0.03)
    r2 = _base_rhs(sc2, F, Tq)
    @test iszero(r2.tree.root_com_bus)
    base_args = engine_args(sc2; mission_time = 10.0)
    p2 = SM.ODEParams(n_sats = 1, args = base_args); SE._initialize_runtime_env_config!(p2); SE._initialize_articulated_runtimes!(p2)
    u2 = SE.build_initial_conditions(base_args); sv = u2.sc[1]; art = p2.shared_buffers.articulated_runtimes[1]
    b2 = AB.ArticulatedBaseState(collect(sv.pos), collect(sv.vel), collect(sv.q), collect(sv.ω))
    _, α2, _ = AB.articulated_dynamics!(art.ws, art.tree, b2, sv.joint_q, sv.joint_qd, F, Tq, r -> zeros(3); root_mass = sv.mass - art.moving_mass)
    @test collect(r2.du.ω) == collect(α2)
end

# ---------------------------------------------------------------------------
# allocations of the articulated RHS
# ---------------------------------------------------------------------------

function _art_rhs_alloc(art, p, u, du, forces, torques, effectors)
    return @allocated SE._assign_articulated_rhs!(du.sc[1], u.sc[1], art, p, 1, 0.0, forces, torques, 0.0, effectors)
end

@testset "articulated RHS is allocation-free" begin
    for eff in (SM.InverseSquaredGravityModel(), SM.InverseSquaredJ2GravityModel())
        args = engine_args(chain_spacecraft(); mission_time = 10.0, effectors = (eff,))
        p = SM.ODEParams(n_sats = 1, args = args)
        SE._initialize_runtime_env_config!(p)
        SE._initialize_articulated_runtimes!(p)
        u = SE.build_initial_conditions(args)
        du = zero(u)
        art = p.shared_buffers.articulated_runtimes[1]
        forces = MVector{3, Float64}(0.1, 0.2, 0.3); torques = MVector{3, Float64}(0.01, 0.0, 0.02)
        SE._assign_articulated_rhs!(du.sc[1], u.sc[1], art, p, 1, 0.0, forces, torques, 0.0, (eff,))    # warm-up
        alloc = _art_rhs_alloc(art, p, u, du, forces, torques, (eff,))
        @info "articulated RHS allocation" effector = nameof(typeof(eff)) alloc
        @test alloc == 0
        # gravity helper agrees with the effector acceleration at the root
        g = SE._gravity_only_acceleration_ii(p, SVector(7.0e6, 1.0e5, -2.0e5), 0.0)
        @test norm(g) > 1.0 && g[1] < 0
    end
end

# ---------------------------------------------------------------------------
# (g) checkpoint and resume
# ---------------------------------------------------------------------------

@testset "checkpoint resume continues an articulated run" begin
    sc() = chain_spacecraft(ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(30.0, -20.0, 10.0); q = rotq([1, -1, 2], 0.9), ang_vel = SVector(0.05, 0.1, -0.04)))
    function with_settings(a; overrides...)
        s = a.simulation_settings
        names = fieldnames(typeof(s))
        vals = NamedTuple{names}(map(n -> getfield(s, n), names))
        return SM.SimConfig._with_configuration(a; simulation_settings = SM.SimulationSettings(; merge(vals, overrides)...))
    end
    dir = mktempdir()
    cfg(mission; kw...) = with_settings(engine_args(sc(); mission_time = mission, data_rate = 5.0, solver_mode = :dp8, tol = tight_tolerances(dt_max = 0.1), results_directory = dir); kw...)
    first_leg = SpaceAGORA.run_simulation(cfg(40.0; checkpoint_enabled = true, checkpoint_interval_s = 20.0); return_results = true)
    resumed = SpaceAGORA.run_simulation(cfg(80.0; checkpoint_enabled = true, checkpoint_interval_s = 20.0, resume_from_checkpoint = true); return_results = true)
    straight = SpaceAGORA.run_simulation(cfg(80.0); return_results = true)
    @test resumed.table.time[end] ≈ 80.0
    cols = ["sc1_joint_q_$k" for k in 1:6]
    cols = vcat(cols, ["sc1_joint_qd_$k" for k in 1:5], ["sc1_pos_$k" for k in 1:3], ["sc1_q_$k" for k in 1:4])
    diff = maximum(abs(resumed.table[end, c] - straight.table[end, c]) / max(1.0, abs(straight.table[end, c])) for c in cols)
    @info "resume vs straight, max relative final-state difference" diff
    @test diff < 1e-8
    @test resumed.table.time[1] >= 40.0 - 1e-9          # the resumed run did not restart from t = 0
end

end # module
