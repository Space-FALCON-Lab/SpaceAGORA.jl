module LivePoseLoadsTests

# Articulated spacecraft with SimulationSettings.articulated_live_pose_loads: per-link aero and facet SRP at
# the live link poses, per-body gravity gradient, the setup refusals, and bit-identity with the switch off.
# Analytic and limit checks; no GRAM or SPICE assets.

using Test
using LinearAlgebra
using StaticArrays
using DataFrames
using SpaceAGORA
using SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const AB = SM.ArticulatedBody
const SE = SpaceAGORA.SimulationEngine

# Allocation counts are only meaningful without coverage instrumentation.
const _ALLOC_CHECKS = Base.JLOptions().code_coverage == 0

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
const PLANET = SM.make_no_gram_planet(:earth)
const MU = Float64(PLANET.μ)

function add_facet!(link; area, normal_vector, cp, ρ, δ)
    SM.add_facet!(link, SM.Facet(area=area, normal_vector=MVector{3, Float64}(normal_vector...), cp=MVector{3, Float64}(cp...), ρ=ρ, δ=δ))
    return link
end

# Dense atmosphere at the test altitude so the aerodynamic loads are well above roundoff.
const RHO_REF = 1.0e-8
const DENSITY = SM.ExponentialAtmosphereModel(RHO_REF, 6.2e5, 6.0e4; temperature_k=1000.0, valid_min_altitude_m=0.0, valid_max_altitude_m=2.0e6)

function tolerances(; dt_max=0.05, rel=1e-12)
    return SM.IntegrationTolerances(
        reltol_orbit=rel, abstol_orbit=1e-9, reltol_atmosphere=rel, abstol_atmosphere=1e-9,
        reltol_quaternion=rel, abstol_quaternion=1e-12, reltol_mass=rel, abstol_mass=1e-12,
        reltol_angular_rate=rel, abstol_angular_rate=1e-12, dt_max_orbit=dt_max, dt_max_atmosphere=dt_max,
    )
end

function with_live(args, flag::Bool)
    st = args.simulation_settings
    names = fieldnames(typeof(st))
    vals = NamedTuple{names}(map(n -> getfield(st, n), names))
    return SM.SimConfig._with_configuration(args; simulation_settings=SM.SimulationSettings(; merge(vals, (articulated_live_pose_loads=flag,))...))
end

function engine_args(sc; mission_time, data_rate=0.05, effectors=(), control=(), live=false, solver_mode=:tsit5,
        tol=tolerances(), density=DENSITY, results_directory=mktempdir())
    base = make_example_config(
        planet=PLANET, spacecraft=sc, mission_time=mission_time, initial_time=T0,
        dynamic_effectors=effectors, density_model=density, ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=true, keplerian=true, verbose=false, results=false, results_directory=results_directory,
        solver_config=SM.SolverConfig(solver_mode=solver_mode),
    )
    args = SM.SimConfig._with_configuration(base;
        integration_tolerances=tol,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=true, num_steps_to_save=1000, data_rate=data_rate),
        control_model=SM.ControlModel(control_effectors=control, control_rates=Float64[]),
    )
    return with_live(args, live)
end

# One articulated spacecraft's RHS machinery, set up the way the engine does it.
function rhs_setup(sc, effectors; live, density=DENSITY, mission_time=10.0)
    args = engine_args(sc; mission_time=mission_time, effectors=effectors, live=live, density=density)
    p = SM.ODEParams(n_sats=1, args=args)
    SE._initialize_runtime_env_config!(p)
    SE._initialize_articulated_runtimes!(p)
    u = SE.build_initial_conditions(args)
    return (p=p, u=u, du=zero(u), art=p.shared_buffers.articulated_runtimes[1], args=args)
end

# Fixed Sun for the engine path: the facet SRP sample reads a position table (the same cache the engine fills
# from SPICE), so no kernels are needed. The table brackets any epoch the tests evaluate.
function with_sun!(S, sun)
    et = S.p.shared_buffers.et_start[]
    S.p.shared_buffers.srp_sun_ephemeris_cache[] = SM.SRPSunEphemerisCache([et - 1.0e6, et + 1.0e6], SVector{3, Float64}[sun, sun])
    return S
end

function eval_rhs!(S, effectors; forces=(0.0, 0.0, 0.0), torques=(0.0, 0.0, 0.0), t=0.0)
    SE._assign_articulated_rhs!(S.du.sc[1], S.u.sc[1], S.art, S.p, 1, t,
        MVector{3, Float64}(forces), MVector{3, Float64}(torques), 0.0, effectors)
    return S
end

# IC with the velocity along +z: the planet-relative flow is mostly along +z (a small -y part from the
# rotating atmosphere), so a hinge about y turns a panel through the flow angle one-for-one.
const IC_Z = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 0.0, 7.5e3);
    q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, 0.0))

# Planet-relative flow velocity at the root COM in inertial axes (the atmosphere is calm).
function flow_ii(S, t=0.0)
    sv = S.u.sc[1]
    pf = SE.sample_planet_frame_with_lpi((pos_ii=SVector{3, Float64}(sv.pos), vel_ii=SVector{3, Float64}(sv.vel)),
        PLANET, SE._planet_lpi_at_engine(S.p, t))
    return pf.l_pi' * pf.vel_pp
end

function rho_at_root(S, t=0.0)
    sv = S.u.sc[1]
    x = SM.StateSample(SVector{3, Float64}(sv.pos), SVector{3, Float64}(sv.vel), 1.0)
    return SE.sample_atmosphere(x, S.p, 1, t; write_buffers=false).rho_kg_m3
end

# Constant-model drag on a link of area `A` and attitude `R_l` (link -> inertial) in the flow `v`.
function expected_constant_drag(ρ, v, R_l, A)
    vb = R_l' * v
    α = atan(vb[1], vb[3])
    a = abs(α)
    CD = 2 * (2.2 - 0.8) / pi * min(a, pi - a) + 0.8
    return 0.5 * ρ * dot(v, v) * CD * A * (-v / norm(v)), CD
end

# Two identical panels at +-x of the bus, hinges about y, panel k deflected by `θ` (rad).
function two_panel_sc(; θR=0.0, θL=0.0, k=0.0, c=0.0, bus_m=20.0, panel_m=2.0, ic=IC_Z, cop=(0.0, 0.0, 0.0), prop=0.0,
        panel_area=1.0, bus_area=1.0, axisL=[0, 1, 0])
    bus = mklink(root=true, m=bus_m, dims=(1.0, 1.0, 1.0), ref_area=bus_area)
    L = mklink(m=panel_m, dims=(1.0, 0.05, 0.5), r=(-1.1, 0.0, 0.0), ref_area=panel_area, cop_offset_b=MVector{3, Float64}(cop...))
    R = mklink(m=panel_m, dims=(1.0, 0.05, 0.5), r=(1.1, 0.0, 0.0), ref_area=panel_area, cop_offset_b=MVector{3, Float64}(cop...))
    jL = mkjoint(bus, L, (-0.5, 0.0, 0.0); joint_type=:hinge, axis=axisL, stiffness=k, damping=c, initial_q=θL)
    jR = mkjoint(bus, R, (0.5, 0.0, 0.0); joint_type=:hinge, axis=[0, 1, 0], stiffness=k, damping=c, initial_q=θR)
    return SM.SpacecraftModel(; joints=[jL, jR], links=[bus, L, R], root=bus, initial_condition=ic, prop_mass=prop)
end

# System net torque about the root COM of the per-body wrenches in the workspace.
function net_torque_about_root(ws, nb)
    τ = zero(SVector{3, Float64})
    for b in 1:nb
        τ += cross(ws.kpos[b] - ws.kpos[1], ws.ext_force[b]) + ws.ext_torque[b]
    end
    return τ
end

# ---------------------------------------------------------------------------
# Switch and default
# ---------------------------------------------------------------------------

@testset "the switch defaults to off and reaches the runtime flag" begin
    @test SM.SimulationSettings().articulated_live_pose_loads === false
    sc = two_panel_sc()
    off = rhs_setup(sc, (SM.AerodynamicCoefficientConstant(),); live=false)
    on = rhs_setup(sc, (SM.AerodynamicCoefficientConstant(),); live=true)
    @test !off.p.shared_buffers.articulated_live_loads[]
    @test on.p.shared_buffers.articulated_live_loads[]
    # a rigid spacecraft never turns the flag on
    bus_r = mklink(root=true, m=10.0, dims=(1.0, 1.0, 1.0))
    rigid = SM.SpacecraftModel(; links=[bus_r], root=bus_r, initial_condition=IC_Z)
    args = engine_args(rigid; mission_time=10.0, live=true)
    p = SM.ODEParams(n_sats=1, args=args)
    SE._initialize_articulated_runtimes!(p)
    @test !p.shared_buffers.articulated_live_loads[]
    # the settings copy helpers carry the field
    @test SpaceAGORA.with_visualization_scene(on.args, true).simulation_settings.articulated_live_pose_loads
end

# ---------------------------------------------------------------------------
# Bit-identity with the switch off
# ---------------------------------------------------------------------------

# Bit patterns of du for identity_rhs(), recorded from the code before the live-pose loads existed (PR #235 tip):
# two hinged panels, point-mass gravity, explicit base loads. Only IEEE-exact operations (no libm) are involved.
const RHS_235 = UInt64[0x403e000000000000, 0x40bce80000000000, 0x4024000000000000, 0xc02045196ca772cd, 0xbfbbaa04a9c56332, 0x3fbad146452c575a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x3f40624dd2f1a9fc, 0xbf50624dd2f1a9fc, 0x3f589374bc6a7efa, 0x8000000000000000, 0x3f5707e953e814a5, 0x3fcc219fca244dbc, 0x3f63beac1e780a9f, 0x0000000000000000, 0x0000000000000000, 0x3fec0a2d54c72be8, 0xc0011115041d63f0]

function bits(v)
    return [reinterpret(UInt64, x) for x in Vector{Float64}(v)]
end

function identity_rhs(; live)
    sc = two_panel_sc(θR=0.3, θL=-0.2, k=5.0, c=0.1, ic=SM.CartesianInitialCondition(SVector(7.0e6, 1.0e5, -2.0e5), SVector(30.0, 7.4e3, 10.0);
        q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(1e-3, -2e-3, 3e-3)), prop=1.0)
    S = rhs_setup(sc, (SM.InverseSquaredGravityModel(),); live=live)
    eval_rhs!(S, (SM.InverseSquaredGravityModel(),); forces=(0.1, 0.2, 0.3), torques=(0.01, 0.0, 0.02))
    return S
end

@testset "switch off: the live-pose code is not entered" begin
    a = identity_rhs(live=false)
    b = identity_rhs(live=false)
    @test bits(a.du.sc[1]) == bits(b.du.sc[1])
    @test bits(a.du.sc[1]) == RHS_235
    # With the switch on but no link-capable effector and no gravity-gradient request, the live path adds nothing:
    # zero external wrenches, so the RHS matches the off path to roundoff (not bits: the on path takes the
    # explicit per-body wrench route).
    c = identity_rhs(live=true)
    @test all(iszero, c.art.ws.ext_force) && all(iszero, c.art.ws.ext_torque)
    @test maximum(abs, collect(c.du.sc[1]) .- collect(a.du.sc[1])) <= 1e-12 * maximum(abs, collect(a.du.sc[1]))
end

# ---------------------------------------------------------------------------
# Per-link drag: symmetric panels, deflected panel
# ---------------------------------------------------------------------------

@testset "symmetric two-panel drag gives zero net torque, equal forces" begin
    eff = (SM.AerodynamicCoefficientConstant(),)
    S = rhs_setup(two_panel_sc(), eff; live=true)
    eval_rhs!(S, eff)
    ws = S.art.ws
    v = flow_ii(S)
    ρ = rho_at_root(S)
    F_exp, CD = expected_constant_drag(ρ, v, Matrix(1.0I, 3, 3), 1.0)
    @test CD ≈ 0.8 rtol = 1e-2                 # near-axial flow: the constant model's CD at zero incidence
    tree = S.art.tree
    @test tree.nb == 3
    # panels (bodies 2, 3): force on each equals the analytic per-link drag (planet rotation at the link COM
    # changes the airspeed by Ω|d|/v ~ 1e-8, hence the tolerance), and no torque (COP at the link COM, no lever)
    for b in 2:3
        @test ws.ext_force[b] ≈ F_exp rtol = 1e-6
        @test norm(ws.ext_torque[b]) <= 1e-6 * norm(F_exp) * 1e-3
    end
    @test norm(ws.ext_force[2] - ws.ext_force[3]) <= 1e-6 * norm(F_exp)
    # root: its own bus drag only
    @test ws.ext_force[1] ≈ F_exp rtol = 1e-6
    # net torque about the root COM: levers +-1.1 m cancel
    τnet = net_torque_about_root(ws, tree.nb)
    @info "symmetric panels: net torque about the root COM" τnet scale = 1.1 * norm(F_exp)
    @test norm(τnet) <= 1e-6 * 1.1 * norm(F_exp)
    # the joint generalized forces on the mirrored panels (hinge axes both +y, levers +-x) are equal and opposite
    qdd = S.du.sc[1].joint_qd
    @test abs(qdd[1]) > 0
    @test abs(qdd[1] + qdd[2]) <= 1e-6 * abs(qdd[1])
end

@testset "deflected panel: analytic force and torque, joint generalized force" begin
    θ0 = 0.4
    cop = (0.0, 0.0, 0.2)
    eff = (SM.AerodynamicCoefficientConstant(),)
    S = rhs_setup(two_panel_sc(θR=θ0, cop=cop), eff; live=true)
    eval_rhs!(S, eff)
    ws = S.art.ws
    v = flow_ii(S)
    ρ = rho_at_root(S)
    Rl = Matrix(AB._rotmat(rotq([0, 1, 0], θ0)))          # panel R: link frame -> inertial
    F_R, CD_R = expected_constant_drag(ρ, v, Rl, 1.0)
    @test CD_R ≈ 0.8 + 2 * 1.4 / pi * θ0 rtol = 1e-3       # alpha = -theta0 (flow along +z), folded
    @test ws.ext_force[3] ≈ F_R rtol = 1e-6
    τ_exp = cross(Rl * collect(cop), F_R)                  # about the panel COM (= its body COM): lever d = 0
    @test ws.ext_torque[3] ≈ τ_exp rtol = 1e-6
    # the undeflected panel L has CD(0) = 0.8-ish and its own COP torque
    F_L, _ = expected_constant_drag(ρ, v, Matrix(1.0I, 3, 3), 1.0)
    @test ws.ext_force[2] ≈ F_L rtol = 1e-6
    @test norm(F_R) > norm(F_L) * 1.3                     # the deflected panel really carries more drag
    # hinge generalized force on R about its hinge: axis . (tau + (x_panel - x_hinge) x F)
    hinge = ws.kpos[1] + SVector(0.5, 0.0, 0.0)
    Q = dot([0.0, 1.0, 0.0], ws.ext_torque[3] + cross(ws.kpos[3] - hinge, ws.ext_force[3]))
    @test isfinite(Q)
    # the root-body drag does not depend on the panels: bus only (area 1, CD 0.8)
    @test ws.ext_force[1] ≈ F_L rtol = 1e-6
end

@testset "a fixed-merged link carries its wrench to the parent body, about that body's COM" begin
    θ0 = 0.3
    bus = mklink(root=true, m=20.0, dims=(1.0, 1.0, 1.0), ref_area=0.0)
    P = mklink(m=2.0, dims=(1.0, 0.05, 0.5), r=(1.1, 0.0, 0.0), ref_area=0.0)
    E = mklink(m=1.0, dims=(0.3, 0.3, 0.3), r=(1.9, 0.0, 0.4), ref_area=0.7)
    jh = mkjoint(bus, P, (0.5, 0.0, 0.0); joint_type=:hinge, axis=[0, 1, 0], stiffness=0.0, damping=0.0, initial_q=θ0)
    jf = mkjoint(P, E, (1.6, 0.0, 0.2))
    sc = SM.SpacecraftModel(; joints=[jh, jf], links=[bus, P, E], root=bus, initial_condition=IC_Z)
    eff = (SM.AerodynamicCoefficientConstant(),)
    S = rhs_setup(sc, eff; live=true)
    eval_rhs!(S, eff)
    tree = S.art.tree
    ws = S.art.ws
    @test tree.nb == 2 && tree.body_of_link == [1, 2, 2]
    v = flow_ii(S)
    ρ = rho_at_root(S)
    R_b = Matrix(AB._rotmat(rotq([0, 1, 0], θ0)))
    F_E, _ = expected_constant_drag(ρ, v, R_b, 0.7)        # E has the identity link attitude in the panel frame
    com_panel = (2.0 * [1.1, 0, 0] + 1.0 * [1.9, 0, 0.4]) / 3.0
    c_E = [1.9, 0.0, 0.4] - com_panel                      # E relative to the panel group's composite COM, panel frame
    @test collect(tree.link_com_in_body[3]) ≈ c_E atol = 1e-12
    @test ws.ext_force[2] ≈ F_E rtol = 1e-6
    @test ws.ext_torque[2] ≈ cross(R_b * c_E, F_E) rtol = 1e-6   # lever about the parent body's COM, not the bus origin
    @test norm(ws.ext_force[1]) == 0.0 && norm(ws.ext_torque[1]) == 0.0      # root drag area is zero
end

# ---------------------------------------------------------------------------
# Relative wind: omega x r
# ---------------------------------------------------------------------------

@testset "the link airspeed includes omega x r (aerodynamic damping)" begin
    ω = SVector(0.0, 0.02, 0.0)              # about y: the panels at +-x move along -+z, i.e. along the flow
    ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 0.0, 7.5e3); q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=ω)
    eff = (SM.AerodynamicCoefficientConstant(),)
    S = rhs_setup(two_panel_sc(ic=ic), eff; live=true)
    eval_rhs!(S, eff)
    ws = S.art.ws
    # velocity of the right panel COM: omega x (1.1, 0, 0) = (0, 0, -0.022) -> less relative speed, less drag
    Fr, Fl = ws.ext_force[3], ws.ext_force[2]
    @test norm(Fr) < norm(Fl)
    ρ = rho_at_root(S)
    v = flow_ii(S)
    for (b, sgn) in ((3, -1.0), (2, +1.0))
        vl = v + SVector(0.0, 0.0, sgn * 1.1 * ω[2])
        F, _ = expected_constant_drag(ρ, vl, Matrix(1.0I, 3, 3), 1.0)
        @test ws.ext_force[b] ≈ F rtol = 1e-6
    end
    # restoring: the net torque about the root COM opposes the rotation (omega . tau < 0)
    @test dot(ω, net_torque_about_root(ws, 3)) < 0
end

# ---------------------------------------------------------------------------
# Stiff joints reproduce the rigid run
# ---------------------------------------------------------------------------

@testset "stiff joints with live-pose loads reproduce the rigid run" begin
    ω0 = SVector(0.0, 0.0, 2e-3)
    ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 0.0, 7.5e3); q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=ω0)
    # k sized so the static deflection tau/k is <= 1e-9 rad: tau <= (aero lever) ~ F * 0.5 m ~ 0.05 N m
    k = 1.0e8
    function build(fixed)
        bus = mklink(root=true, m=20.0, dims=(1.0, 1.0, 1.0))
        L = mklink(m=2.0, dims=(1.0, 0.05, 0.5), r=(-1.1, 0.0, 0.0))
        R = mklink(m=2.0, dims=(1.0, 0.05, 0.5), r=(1.1, 0.0, 0.0))
        kw = fixed ? (;) : (; joint_type=:hinge, axis=[0, 1, 0], stiffness=k, damping=2 * sqrt(k * 0.2) * 0.5)
        joints = [mkjoint(bus, L, (-0.5, 0.0, 0.0); kw...), mkjoint(bus, R, (0.5, 0.0, 0.0); kw...)]
        return SM.SpacecraftModel(; joints=joints, links=[bus, L, R], root=bus, initial_condition=ic,
            inertia_tensor=bus.inertia)
    end
    # composite inertia of the all-locked configuration
    tree_locked = SM.build_articulated_tree(build(true))
    Jc = tree_locked.inertia[1]
    function rigid_sc()
        sc = build(true)
        sc.inertia_tensor = Jc
        return sc
    end
    tend = 120.0
    # Gravity gradient in both runs: the rigid run applies it from the composite inertia, the live run per body
    # (per-body point gravity plus the per-body intrinsic term); they agree to O(d/r) ~ 1e-7.
    eff = (SM.InverseSquaredGravityModel(gravity_gradient=true), SM.AerodynamicCoefficientConstant())
    tol = tolerances(dt_max=0.5, rel=1e-12)
    r_rigid = SpaceAGORA.run_simulation(engine_args(rigid_sc(); mission_time=tend, data_rate=tend, effectors=eff, solver_mode=:rodas5p, tol=tol); return_results=true)
    r_live = SpaceAGORA.run_simulation(engine_args(build(false); mission_time=tend, data_rate=tend, effectors=eff, live=true, solver_mode=:rodas5p, tol=tol); return_results=true)
    dpos = norm([r_rigid.table.sc1_pos_1[end] - r_live.table.sc1_pos_1[end], r_rigid.table.sc1_pos_2[end] - r_live.table.sc1_pos_2[end],
        r_rigid.table.sc1_pos_3[end] - r_live.table.sc1_pos_3[end]])
    dvel = norm([r_rigid.table.sc1_vel_1[end] - r_live.table.sc1_vel_1[end], r_rigid.table.sc1_vel_2[end] - r_live.table.sc1_vel_2[end],
        r_rigid.table.sc1_vel_3[end] - r_live.table.sc1_vel_3[end]])
    # Neglected terms (derived): live airspeed adds |omega||r_link|/|v| = 2e-3*1.1/7.5e3 ~ 3e-7 of the drag at the
    # panels (and the planet rotation seen at the link COM, 1e-8); the drag displacement after tend is
    # 0.5 (F/m) t^2 with F ~ 0.5 rho v^2 CD A_total.
    S0 = rhs_setup(rigid_sc(), eff; live=false)
    ρ = rho_at_root(S0)
    v = norm(flow_ii(S0))
    F = 0.5 * ρ * v^2 * 2.2 * 3.0
    drag_disp = 0.5 * F / (20 + 4) * tend^2
    rel = 2e-3 * 1.1 / v + 7.3e-5 * 1.1 / v
    bound_pos = 50 * rel * drag_disp + 1e-6
    @info "stiff joints vs rigid" dpos dvel bound_pos drag_disp
    @test drag_disp > 1e-2                                   # the drag matters in this run
    @test dpos <= bound_pos
    # Attitude: the only torque the rigid run lacks is the aerodynamic damping from omega x r at the panels,
    # tau_d = 2 r^2 omega F_panel / v about z; its attitude effect is dq ~ 0.5 * (0.5 tau_d / Izz t^2) (quaternion
    # half-angle). The static joint deflection (<= 1e-9 rad) is far below it.
    F_panel = 0.5 * ρ * v^2 * 0.8 * 1.0
    τ_d = 2 * 1.1^2 * norm(ω0) * F_panel / v
    dq_est = 0.5 * 0.5 * τ_d / Jc[3, 3] * tend^2
    q_rigid = [r_rigid.table[end, "sc1_q_$i"] for i in 1:4]
    q_live = [r_live.table[end, "sc1_q_$i"] for i in 1:4]
    @info "stiff joints vs rigid, attitude" dq = norm(q_rigid - q_live) dq_est
    @test norm(q_rigid - q_live) <= 2 * dq_est
end

# ---------------------------------------------------------------------------
# Gravity-gradient libration of a hinged dumbbell
# ---------------------------------------------------------------------------

function dumbbell_run(; gg::Bool, mission_time, φ0=0.02)
    r0 = 7.0e6
    n = sqrt(MU / r0^3)
    ic = SM.CartesianInitialCondition(SVector(r0, 0.0, 0.0), SVector(0.0, sqrt(MU / r0), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, n))
    m_link = 1.0
    d = 1.0
    L = 2.4
    bus = mklink(root=true, m=1.0e8, dims=(1.0, 1.0, 1.0))          # cube: spherical inertia
    rod = mklink(m=m_link, dims=(L, 0.01, 0.01), r=(d, 0.0, 0.0))
    j = mkjoint(bus, rod, (0.0, 0.0, 0.0); joint_type=:hinge, axis=[0, 0, 1], stiffness=0.0, damping=0.0, initial_q=φ0)
    sc = SM.SpacecraftModel(; joints=[j], links=[bus, rod], root=bus, initial_condition=ic)
    eff = (SM.InverseSquaredGravityModel(gravity_gradient=gg),)
    args = engine_args(sc; mission_time=mission_time, data_rate=5.0, effectors=eff, live=true, solver_mode=:dp8,
        tol=tolerances(dt_max=5.0, rel=1e-12), density=SM.NoAtmosphereModel())
    res = SpaceAGORA.run_simulation(args; return_results=true)
    Ic = m_link * (L^2 + 0.01^2) / 12                                  # about the hinge axis (z), through the rod COM
    Ia = m_link * (0.01^2 + 0.01^2) / 12                               # along the rod
    return (table=res.table, n=n, d=d, m=m_link, Ic=Ic, Ia=Ia)
end

function libration_frequency(tbl)
    zc = Float64[]
    ts, ys = tbl.time, tbl.sc1_joint_q_1
    for i in 2:length(ts)
        ys[i - 1] * ys[i] < 0 && push!(zc, ts[i - 1] - ys[i - 1] * (ts[i] - ts[i - 1]) / (ys[i] - ys[i - 1]))
    end
    nz = length(zc)
    return nz, pi * (nz - 1) / (zc[end] - zc[1])                       # half period per crossing
end

@testset "hinged dumbbell libration at sqrt(3) n, and it fails without the per-body gravity gradient" begin
    n0 = sqrt(MU / 7.0e6^3)
    tend = 12 * 2π / (sqrt(3) * n0)
    # With the gravity gradient: sqrt(3 (C - A) / C) n with C the pitch inertia about the hinge, A the inertia along the rod.
    on = dumbbell_run(gg=true, mission_time=tend)
    C = on.Ic + on.m * on.d^2
    ω_exact = sqrt(3 * (C - on.Ia) / C) * on.n
    nz, ω_on = libration_frequency(on.table)
    err_on = abs(ω_on - sqrt(3) * on.n) / (sqrt(3) * on.n)
    err_exact = abs(ω_on - ω_exact) / ω_exact
    @info "dumbbell libration with per-body GG" nz ω_on_over_n = ω_on / on.n sqrt3 = sqrt(3) err_vs_sqrt3 = err_on err_vs_exact = err_exact
    @test nz >= 20
    @test err_on < 2e-3
    @test err_exact < 2e-3
    # Without it (gravity_gradient=false), only the point-mass term acts: sqrt(3 m d^2 / (Ic + m d^2)) n, >= 10 percent lower.
    off = dumbbell_run(gg=false, mission_time=tend)
    nz0, ω_off = libration_frequency(off.table)
    ω_missing = sqrt(3 * off.m * off.d^2 / C) * off.n
    @info "dumbbell libration without GG" ω_off_over_n = ω_off / off.n predicted_over_n = ω_missing / off.n shortfall = 1 - ω_off / ω_on
    @test (ω_on - ω_off) / ω_on > 0.10
    @test abs(ω_off - ω_missing) / ω_missing < 5e-3
end

@testset "gravity-gradient request: per-body torque, the configured-inertia versions are dropped" begin
    ic = SM.CartesianInitialCondition(SVector(7.0e6, 1.0e5, -2.0e5), SVector(0.0, 7.5e3, 0.0); q=rotq([1, 2, 3], 0.7), ang_vel=SVector(0.0, 0.0, 0.0))
    sc = two_panel_sc(θR=0.2, θL=-0.1, ic=ic)
    for eff in ((SM.InverseSquaredGravityModel(gravity_gradient=true),),
                (SM.InverseSquaredGravityModel(), SM.GravityGradientTorqueModel()),
                (SM.InverseSquaredJ2GravityModel(gravity_gradient=true),))
        S = rhs_setup(sc, eff; live=true)
        eval_rhs!(S, eff)
        ws = S.art.ws
        tree = S.art.tree
        for b in 1:tree.nb
            R = AB._rotmat(ws.kquat[b])
            x = SVector{3, Float64}(S.u.sc[1].pos) + ws.kpos[b]
            r̂ = x / norm(x)
            τ = 3 * MU / norm(x)^3 * cross(r̂, (R * tree.inertia[b] * R') * r̂)
            @test ws.ext_torque[b] ≈ τ rtol = 1e-12
        end
        # the stand-alone GG effector is not applied on the root as well: it is filtered from the base effectors
        @test SE._live_base_effectors(eff) == ()
    end
    # off: v1 semantics (effector kept on the root; a gravity_gradient flag is dropped with its force)
    @test SE._nongravity_effectors((SM.InverseSquaredGravityModel(), SM.GravityGradientTorqueModel())) isa Tuple{SM.GravityGradientTorqueModel}
    # no request, no per-body torque
    S = rhs_setup(sc, (SM.InverseSquaredGravityModel(),); live=true)
    eval_rhs!(S, (SM.InverseSquaredGravityModel(),))
    @test all(iszero, S.art.ws.ext_torque)
end

@testset "gravity-gradient request on an articulated spacecraft is refused unless the switch is on" begin
    sc = two_panel_sc(θR=0.2, θL=-0.1)
    for eff in ((SM.InverseSquaredGravityModel(gravity_gradient=true),),
                (SM.ConstantGravityModel(gravity_gradient=true),),
                (SM.InverseSquaredJ2GravityModel(gravity_gradient=true),),
                (SM.InverseSquaredGravityModel(), SM.GravityGradientTorqueModel()))
        err = try
            SE._validate_articulated_spacecraft!(engine_args(sc; mission_time=10.0, effectors=eff, live=false), :tsit5)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("articulated_live_pose_loads=true", err.msg)
        @test occursin(string(nameof(typeof(eff[end]))), err.msg)
        # the same configuration is accepted and runs with the switch on
        @test SE._validate_articulated_spacecraft!(engine_args(sc; mission_time=10.0, effectors=eff, live=true), :tsit5) === nothing
    end
    eff = (SM.InverseSquaredGravityModel(gravity_gradient=true),)
    @test_throws ArgumentError SpaceAGORA.run_simulation(engine_args(sc; mission_time=5.0, data_rate=1.0, effectors=eff, live=false); return_results=true)
    res = SpaceAGORA.run_simulation(engine_args(sc; mission_time=5.0, data_rate=1.0, effectors=eff, live=true); return_results=true)
    @test all(isfinite, res.table.sc1_pos_1)
    # no request: accepted with the switch off; a disabled GravityGradientTorqueModel is not a request
    @test SE._validate_articulated_spacecraft!(engine_args(sc; mission_time=10.0, effectors=(SM.InverseSquaredGravityModel(gravity_gradient=false),), live=false), :tsit5) === nothing
    @test SE._validate_articulated_spacecraft!(engine_args(sc; mission_time=10.0, effectors=(SM.GravityGradientTorqueModel(gravity_gradient=false),), live=false), :tsit5) === nothing
    # rigid spacecraft are unaffected
    bus_r = mklink(root=true, m=10.0, dims=(1.0, 1.0, 1.0))
    rigid = SM.SpacecraftModel(; links=[bus_r], root=bus_r, initial_condition=IC_Z, inertia_tensor=bus_r.inertia)
    @test SE._validate_articulated_spacecraft!(engine_args(rigid; mission_time=10.0, effectors=eff, live=false), :tsit5) === nothing
end

# ---------------------------------------------------------------------------
# Facet SRP
# ---------------------------------------------------------------------------

const SUN = SVector(1.495978707e11, 0.0, 0.0)
const SRP = SM.FacetSolarRadiationPressureModel()

function srp_link(; facet_normal=(1.0, 0.0, 0.0), cp=(0.0, 0.0, 0.3), area=0.8, ρ=0.0, δ=0.0)
    link = mklink(m=1.0, dims=(0.05, 1.0, 0.5))
    add_facet!(link; area=area, normal_vector=facet_normal, cp=cp, ρ=ρ, δ=δ)
    return link
end

link_env(pos) = SM.EnvironmentSample(PLANET; solar=SM.SolarEphemerisSample(SUN))
lstate(pos, q) = SM.LinkStateSample(1, 2, pos, SVector(0.0, 0.0, 0.0), q, SVector(0.0, 0.0, 0.0))

@testset "facet SRP at a rotated link: cosine law, edge-on, back face, shadow" begin
    link = srp_link()
    pos = SVector(7.0e6, 0.0, 0.0)
    d = norm(SUN - pos)
    P = SRP.p_srp_1au * (SRP.AU_m / d)^2
    ŝ = (SUN - pos) / d
    for θ0 in (0.0, 0.5, 1.2)
        q = rotq([0, 0, 1], θ0)
        F, τ, _, _, _ = SM.link_wrench(SRP, link, lstate(pos, q), link_env(pos), 0.0, nothing, 0)
        R = Matrix(AB._rotmat(q))
        n_ii = R * [1.0, 0.0, 0.0]
        cosθ = dot(n_ii, ŝ)
        F_exp = -P * 0.8 * cosθ * ŝ                         # absorbing facet: rho = delta = 0
        @test F ≈ F_exp rtol = 1e-12
        @test τ ≈ cross(R * [0.0, 0.0, 0.3], F_exp) rtol = 1e-12
    end
    # edge on (normal perpendicular to the Sun line): force and torque vanish to roundoff
    q_edge = rotq([0, 0, 1], π / 2)
    F, τ, _, _, _ = SM.link_wrench(SRP, link, lstate(pos, q_edge), link_env(pos), 0.0, nothing, 0)
    @test norm(F) <= 1e-12 * P * 0.8
    @test norm(τ) <= 1e-12 * P * 0.8 * 0.3
    # back face: the facet does not see the Sun, exactly zero
    q_back = rotq([0, 0, 1], π)
    F, τ, _, _, _ = SM.link_wrench(SRP, link, lstate(pos, q_back), link_env(pos), 0.0, nothing, 0)
    @test F == zeros(3) && τ == zeros(3)
    # shadow, per link: the same link behind the planet carries no force; in front it does
    q = rotq([0, 0, 1], 0.5)
    Fs, τs, _, _, _ = SM.link_wrench(SRP, link, lstate(SVector(-7.0e6, 0.0, 0.0), q), link_env(nothing), 0.0, nothing, 0)
    @test Fs == zeros(3) && τs == zeros(3)
    Fl, _, _, _, _ = SM.link_wrench(SRP, link, lstate(pos, q), link_env(nothing), 0.0, nothing, 0)
    @test norm(Fl) > 0
    # a link with no facets contributes nothing
    Fn, τn, _, _, _ = SM.link_wrench(SRP, mklink(m=1.0, dims=(0.05, 1.0, 0.5)), lstate(pos, q), link_env(nothing), 0.0, nothing, 0)
    @test Fn == zeros(3) && τn == zeros(3)
end

@testset "facet SRP through the engine: wrench on the carrying body at the live link pose" begin
    θ0 = 0.5
    bus = mklink(root=true, m=20.0, dims=(1.0, 1.0, 1.0))
    P = mklink(m=2.0, dims=(1.0, 0.05, 0.5), r=(1.1, 0.0, 0.0))
    add_facet!(P; area=0.8, normal_vector=(0.0, 0.0, 1.0), cp=(0.0, 0.0, 0.2), ρ=0.0, δ=0.0)
    j = mkjoint(bus, P, (0.5, 0.0, 0.0); joint_type=:hinge, axis=[0, 1, 0], stiffness=0.0, damping=0.0, initial_q=θ0)
    sc = SM.SpacecraftModel(; joints=[j], links=[bus, P], root=bus, initial_condition=IC_Z)
    eff = (SRP,)
    S = with_sun!(rhs_setup(sc, eff; live=true), SUN)
    eval_rhs!(S, eff)
    ws = S.art.ws
    sv = S.u.sc[1]
    # independent link pose from the library
    lp, lq = AB.articulated_link_poses(S.art.tree, [ws.kpos[b] + SVector{3, Float64}(sv.pos) for b in 1:2], ws.kquat)
    sun = SE.sample_solar_ephemeris(nothing, S.p, 1, 0.0).sun_pos_ii
    R = Matrix(AB._rotmat(lq[2]))
    ŝ = (sun - lp[2]) / norm(sun - lp[2])
    n_ii = R * [0.0, 0.0, 1.0]
    P_srp = SRP.p_srp_1au * (SRP.AU_m / norm(sun - lp[2]))^2 * SM.DynamicEffectors.eclipse_area_calc(lp[2], sun, Float64(PLANET.Rp_e))
    cosθ = dot(n_ii, ŝ)
    F_exp = cosθ > 0 ? -P_srp * 0.8 * cosθ * ŝ : zeros(3)
    @info "engine SRP on the hinged panel" cosθ shadow = P_srp / (SRP.p_srp_1au * (SRP.AU_m / norm(sun - lp[2]))^2)
    @test ws.ext_force[2] ≈ F_exp rtol = 1e-9
    @test ws.ext_torque[2] ≈ cross(R * [0.0, 0.0, 0.2], F_exp) rtol = 1e-9
    @test norm(ws.ext_force[1]) == 0.0                    # the bus carries no facets
end

# ---------------------------------------------------------------------------
# Live link sample: panel-angle control on the root body, refusals
# ---------------------------------------------------------------------------

struct PanelCtl <: SM.AbstractForceTorqueModel
    controlled_panel_links::Vector{Int}
end

@testset "root-body fixed links use the live link attitude; moving-link panel control is refused" begin
    bus = mklink(root=true, m=20.0, dims=(1.0, 1.0, 1.0))
    P = mklink(m=2.0, dims=(1.0, 0.05, 0.5), r=(1.1, 0.0, 0.0))
    E = mklink(m=1.0, dims=(0.3, 0.3, 0.3), r=(0.0, 0.9, 0.0), ref_area=0.5)
    jh = mkjoint(bus, P, (0.5, 0.0, 0.0); joint_type=:hinge, axis=[0, 1, 0], stiffness=0.0, damping=0.0, initial_q=0.2)
    jf = mkjoint(bus, E, (0.0, 0.5, 0.0))
    sc = SM.SpacecraftModel(; joints=[jh, jf], links=[bus, P, E], root=bus, initial_condition=IC_Z)
    S = rhs_setup(sc, (SM.AerodynamicCoefficientConstant(),); live=true)
    eval_rhs!(S, (SM.AerodynamicCoefficientConstant(),))
    ws = S.art.ws
    tree = S.art.tree
    base_pos = SVector{3, Float64}(S.u.sc[1].pos)
    xl0, _ = SE._live_link_sample(ws, tree, sc.links[3], 3, 3, base_pos, SVector(0.0, 0.0, 0.0))
    @test xl0.body == 1
    # a controller rewrites the fixed link's attitude: the live sample follows it
    sc.links[3].q .= rotq([0, 1, 0], 0.3)
    xl1, _ = SE._live_link_sample(ws, tree, sc.links[3], 3, 3, base_pos, SVector(0.0, 0.0, 0.0))
    @test xl1.q_ib ≈ AB._qnormalize(AB._qmul(ws.kquat[1], rotq([0, 1, 0], 0.3))) atol = 1e-14
    @test norm(xl1.q_ib - xl0.q_ib) > 1e-3
    # the moving link's attitude is the kinematic one (joint angle 0.2 about y)
    xp, _ = SE._live_link_sample(ws, tree, sc.links[2], 2, 2, base_pos, SVector(0.0, 0.0, 0.0))
    @test xp.q_ib ≈ AB._qnormalize(AB._qmul(ws.kquat[2], tree.link_q_in_body[2])) atol = 1e-14
    # panel-angle control of a link on a moving body stays refused, with the switch on and off
    for flag in (true, false)
        args = engine_args(sc; mission_time=10.0, effectors=(PanelCtl([1]),), live=flag)   # link 1 is the bus; link 2 is the hinged panel
        @test SE._validate_articulated_spacecraft!(args, :tsit5) === nothing
        args2 = engine_args(sc; mission_time=10.0, effectors=(PanelCtl([2]),), live=flag)
        @test_throws ArgumentError SE._validate_articulated_spacecraft!(args2, :tsit5)
    end
end

@testset "a whole-vehicle (root_only) mesh surrogate is refused with joints moving and the switch on" begin
    coeff = reshape(sqrt(4pi) .* [-1.0, 0.2, 0.1, 0.05, 0.1, -0.07], 6, 1)
    sur = SM.MeshAeroSurrogate(0, 0, coeff, zeros(6, 1), 1.0, 1.0, SVector(0.0, 0.0, 0.0), 1.0, 1.0, 3.0, 20.0)
    whole = SM.AerodynamicCoefficientMeshSurrogate(Dict(1 => sur), true)
    sc = two_panel_sc(θR=0.1)
    @test_throws ArgumentError SE._validate_articulated_spacecraft!(engine_args(sc; mission_time=10.0, effectors=(whole,), live=true), :tsit5)
    @test SE._validate_articulated_spacecraft!(engine_args(sc; mission_time=10.0, effectors=(whole,), live=false), :tsit5) === nothing
end

# ---------------------------------------------------------------------------
# Allocations
# ---------------------------------------------------------------------------

@testset "live-pose RHS is allocation-free with point-mass gravity" begin
    bus = mklink(root=true, m=20.0, dims=(1.0, 1.0, 1.0))
    P = mklink(m=2.0, dims=(1.0, 0.05, 0.5), r=(1.1, 0.0, 0.0))
    add_facet!(P; area=0.8, normal_vector=(0.0, 0.0, 1.0), cp=(0.0, 0.0, 0.2), ρ=0.1, δ=0.3)
    j = mkjoint(bus, P, (0.5, 0.0, 0.0); joint_type=:hinge, axis=[0, 1, 0], stiffness=1.0, damping=0.1, initial_q=0.3)
    ic = SM.CartesianInitialCondition(SVector(7.0e6, 1.0e5, -2.0e5), SVector(10.0, 7.4e3, 20.0); q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(1e-3, 2e-3, -1e-3))
    sc = SM.SpacecraftModel(; joints=[j], links=[bus, P], root=bus, initial_condition=ic)
    for aero in (SM.AerodynamicCoefficientfM(), SM.AerodynamicCoefficientConstant())
        eff = (SM.InverseSquaredGravityModel(gravity_gradient=true), aero, SRP)
        S = with_sun!(rhs_setup(sc, eff; live=true), SUN)
        eval_rhs!(S, eff; forces=(0.1, 0.2, 0.3), torques=(0.01, 0.0, 0.02))      # warm-up
        eval_rhs!(S, eff; forces=(0.1, 0.2, 0.3), torques=(0.01, 0.0, 0.02))
        f = MVector{3, Float64}(0.1, 0.2, 0.3); τ = MVector{3, Float64}(0.01, 0.0, 0.02)
        alloc = @allocated SE._assign_articulated_rhs!(S.du.sc[1], S.u.sc[1], S.art, S.p, 1, 0.0, f, τ, 0.0, eff)
        @info "live-pose RHS allocation" aero = nameof(typeof(aero)) alloc
        @test alloc == 0 skip = !_ALLOC_CHECKS
        @test any(!iszero, S.art.ws.ext_force[2])                                  # the loads are really evaluated
    end
end

end # module
