# Stage 2: a ProximityScene inside `run_simulation` (shadow entries). Included after common.jl.
using Test

# --- configuration helpers ----------------------------------------------------------------------------

default_tol() = SM.IntegrationTolerances(reltol_orbit=1e-12, abstol_orbit=1e-9, reltol_atmosphere=1e-12, abstol_atmosphere=1e-9,
    reltol_quaternion=1e-12, abstol_quaternion=1e-9, reltol_mass=1e-12, abstol_mass=1e-9,
    reltol_angular_rate=1e-12, abstol_angular_rate=1e-9, dt_max_orbit=0.5, dt_max_atmosphere=0.5)

# A configuration over `scs` (two spacecraft by default) with the scene owning `owned` (pairs index => body).
function engine_config(effectors, tend, data_rate; scs=nothing, scene=nothing, owned=(1 => "chaser", 2 => "target"),
        solver=:dp8, orient=false, tol=default_tol(), extra=(;))
    scs === nothing && (scs = SM.SpacecraftModel[rigid_sc(20.0, (1.0, 0.8, 0.6), R1, V1), rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)])
    base = make_example_config(planet=EARTH, spacecraft=scs[1], mission_time=tend, initial_time=T0,
        dynamic_effectors=effectors, density_model=SM.NoAtmosphereModel(), ephemerides_model=EPH,
        orientation_sim=orient, keplerian=true, verbose=false, results=false, results_directory=mktempdir(),
        solver_config=SM.SolverConfig(solver_mode=solver))
    props = scene === nothing ? () : (ProximitySceneDynamics(scene, owned...),)
    return SM.SimConfig._with_configuration(base;
        integration_tolerances=tol,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=tend, orientation_sim=orient, num_steps_to_save=100000, data_rate=data_rate),
        dynamics_model=SM.DynamicsModel(scs, effectors), external_propagators=props, extra...)
end
run_table(cfg; kwargs...) = SpaceAGORA.run_simulation(cfg; return_results=true, kwargs...).table
# A do-nothing periodic callback: its tstops make the solver land exactly on every save time, so saved values are
# the synced shadow states and not dense-output interpolants.
tick_every(period) = (SM.SimulationCallbacks.PeriodicCallback(_ -> nothing, period),)

@testset "scene inside run_simulation: bit-identical to the standalone runner, close to ordinary spacecraft" begin
    effectors = (SM.InverseSquaredJ2GravityModel(),)
    dt = 0.0125
    tend = 600.0; every = 50.0
    scene = make_scene(effectors; dt)
    tbl = run_table(engine_config(effectors, tend, every; scene); extra_callbacks=tick_every(every))
    got = samples_from_table(tbl, every, tend)
    ref_scene = scene_samples(effectors, dt, tend, every)
    @test length(got) == length(ref_scene)
    # The scene steps do not depend on the engine, so the saved shadow entries (a sync lands on every save time)
    # equal the standalone runner's states exactly.
    maxdiff = maximum(maximum(norm.(collect(g) .- collect(r))) for (g, r) in zip(got, ref_scene))
    @info "engine vs standalone runner: max component-wise difference" maxdiff
    @test got == ref_scene
    # And they track the same two bodies flown as ordinary spacecraft within the Stage 1 convergence bounds.
    ordinary = samples_from_table(run_table(engine_config(effectors, tend, every); extra_callbacks=tick_every(every)), every, tend)
    d = max_differences(got, ordinary)
    @info "engine scene vs ordinary spacecraft over $(tend) s" d
    @test d.dr < TOL_POSITION_M && d.dv < TOL_VELOCITY_MPS && d.drel < TOL_RELATIVE_POSITION_M
end

const SE = SpaceAGORA.SimulationEngine
const MU_EARTH = MU

@testset "shadow entries between syncs follow the chief acceleration within a bounded error" begin
    # A coarse scene step (dt = 0.25 s) makes the solver sync between ticks. Saved rows at every tick are then
    # dense-output values of shadow entries extrapolated from a scene state up to one step old. Against the
    # scene's own state at that tick the error is the integral of the difference between the body's and the
    # chief's acceleration (the tidal term, at most 3 n^2 |rho|) over at most one scene step.
    effectors = (SM.InverseSquaredGravityModel(),)
    dt = 0.25; tend = 100.0
    scene = make_scene(effectors; dt)
    got = samples_from_table(run_table(engine_config(effectors, tend, dt; scene)), dt, tend)
    truth = scene_samples(effectors, dt, tend, dt)
    d = max_differences(got, truth)
    n2 = MU_EARTH / norm(R1)^3
    rho = norm(R2 - R1) + 1.0
    # Each of the two contributions (scene lag extrapolated at the sync, then interpolation) is at most
    # (1/2) (3 n^2 rho) dt^2 in position and (3 n^2 rho) dt in velocity; 1.5 is the margin on that sum.
    bound_r = 1.5 * 3 * n2 * rho * dt^2
    bound_v = 1.5 * 3 * n2 * rho * 2dt
    @info "between-sync shadow error" max_position_m = d.dr max_velocity_mps = d.dv bound_position_m = bound_r bound_velocity_mps = bound_v
    @test d.dr < bound_r
    @test d.dv < bound_v
    @test d.dr > 0                         # the rows really are off-tick interpolants, not the exact states
end

# A navigation effector that records the shadow entry of spacecraft 1 whenever it runs.
mutable struct ShadowProbe
    t::Vector{Float64}
    pos::Vector{SVector{3, Float64}}
end
function SM.NavigationHooks.calcNavigationEffect!(m::ShadowProbe, u, p, t::Float64, sat_idx::Int)
    if sat_idx == 1
        push!(m.t, t); push!(m.pos, SVector{3, Float64}(u.sc[1].pos))
    end
    return nothing
end

@testset "guidance and navigation read synced shadow state at GNC ticks" begin
    effectors = (SM.InverseSquaredGravityModel(),)
    dt = 0.05; period = 0.5; tend = 20.0
    probe = ShadowProbe(Float64[], SVector{3, Float64}[])
    scene = make_scene(effectors; dt)
    cfg = engine_config(effectors, tend, 5.0; scene, extra=(; navigation_model=SM.NavigationModel(navigation_effectors=(probe,), navigation_rates=[period])))
    SpaceAGORA.run_simulation(cfg; isolate_state=false)
    @test length(probe.t) >= round(Int, tend / period) - 1
    ref = make_scene(effectors; dt)
    steps_per = round(Int, period / dt)
    ok = true
    for (k, (t, pos)) in enumerate(zip(probe.t, probe.pos))
        k == 1 && continue                     # the initialization call at t = 0 runs before any sync
        n = round(Int, t / dt)
        while ref.n < n; scene_step!(ref); end
        ok &= pos == scene_body_state(ref, 1).r
    end
    @test ok                                   # exactly the scene state at the tick: the sync callback runs first
    # a rate that is not a multiple of the scene step is refused
    bad = engine_config(effectors, tend, 5.0; scene, extra=(; navigation_model=SM.NavigationModel(navigation_effectors=(probe,), navigation_rates=[0.07])))
    @test_throws ArgumentError SE._validate_external_propagators!(bad, :dp8)
    good = engine_config(effectors, tend, 5.0; scene, extra=(; navigation_model=SM.NavigationModel(navigation_effectors=(probe,), navigation_rates=[0.1])))
    @test SE._validate_external_propagators!(good, :dp8) === nothing
end

# --- attitude convention ------------------------------------------------------------------------------

const QMATH = SM.QuaternionMath
quat_angle(a, b) = (c = abs(dot(a, b)); 2 * acos(min(1.0, c)))     # angle between two attitude quaternions

@testset "attitude convention: SpaceAGORA scalar-last q <-> MuJoCo xquat" begin
    using Random
    rng = MersenneTwister(20261007)
    # The mapping is a reordering of components: exact, and its own inverse.
    for _ in 1:100
        q = SVector{4, Float64}(randn(rng, 4))
        @test mujoco_to_sa_quaternion(sa_to_mujoco_quaternion(q)) === q
        @test sa_to_mujoco_quaternion(mujoco_to_sa_quaternion(q)) === q
    end
    @test sa_to_mujoco_quaternion(SVector(0.1, 0.2, 0.3, 0.9)) == SVector(0.9, 0.1, 0.2, 0.3)

    # MuJoCo's own body-to-world matrix (ximat) equals the transpose of SpaceAGORA's rot(q), which is
    # inertial-to-body; the body rate in the free joint is the body-frame rate, and its world-axes image is R*omega.
    ang = 0.0; bad_if_transposed = 0.0
    for _ in 1:20
        q = SVector{4, Float64}(normalize(randn(rng, 4)))
        ω = SVector{3, Float64}(0.1 .* randn(rng, 3))
        st = [SceneBodyState("chaser", R1, V1; q=sa_to_mujoco_quaternion(q), ω=ω), SceneBodyState("target", R2, V2)]
        sc = ProximityScene(; mjcf_xml=TWO_BODIES, dt=0.05, planet=EARTH, gravity_effectors=(SM.InverseSquaredGravityModel(),), initial_states=st)
        Binding.forward!(sc.model, sc.data)
        X = Binding.ximat(sc.model, sc.data)[:, 2]
        R_mj = SMatrix{3, 3, Float64}(X[1], X[4], X[7], X[2], X[5], X[8], X[3], X[6], X[9])    # row-major to matrix
        R_sa = QMATH.rot(q)
        ang = max(ang, maximum(abs.(R_mj - R_sa')))
        bad_if_transposed = max(bad_if_transposed, maximum(abs.(R_mj - R_sa)))
        s = scene_body_state(sc, "chaser")
        @test maximum(abs.(s.ω_world - R_sa' * s.ω)) < 1e-15
    end
    @info "ximat vs rot(q)'" max_abs_difference = ang max_abs_difference_if_not_transposed = bad_if_transposed
    @test ang < 1e-15
    @test bad_if_transposed > 0.1                  # the transpose matters: this check can fail

    # Round trip through the scene (reset, forward kinematics, read back). Quaternions that are exactly unit
    # norm in floating point survive MuJoCo's renormalization bit for bit; the body-frame rate always does.
    exact_unit = (SVector(0.5, 0.5, 0.5, 0.5), SVector(0.5, -0.5, 0.5, 0.5), SVector(0.0, 0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0, 0.0),
        SVector(0.0, 0.6, 0.0, 0.8), SVector(0.6, 0.0, 0.8, 0.0))
    for q in exact_unit
        @test norm(q) == 1.0
        ω = SVector{3, Float64}(randn(rng, 3))
        ic = SM.CartesianInitialCondition(R1, V1; q=q, ang_vel=ω)
        st = [body_state_from_initial_condition("chaser", ic), SceneBodyState("target", R2, V2)]
        sc = ProximityScene(; mjcf_xml=TWO_BODIES, dt=0.05, planet=EARTH, gravity_effectors=(SM.InverseSquaredGravityModel(),), initial_states=st)
        s = scene_body_state(sc, "chaser")
        @test mujoco_to_sa_quaternion(s.q) === q
        @test s.ω === ω
    end
    # Arbitrary unit quaternions come back to a few ulp (MuJoCo renormalizes).
    worst = 0.0
    for _ in 1:50
        q = SM.project_unit_quaternion(SVector{4, Float64}(randn(rng, 4)))
        ω = SVector{3, Float64}(randn(rng, 3))
        st = [SceneBodyState("chaser", R1, V1; q=sa_to_mujoco_quaternion(q), ω=ω), SceneBodyState("target", R2, V2)]
        sc = ProximityScene(; mjcf_xml=TWO_BODIES, dt=0.05, planet=EARTH, gravity_effectors=(SM.InverseSquaredGravityModel(),), initial_states=st)
        s = scene_body_state(sc, "chaser")
        worst = max(worst, maximum(abs.(mujoco_to_sa_quaternion(s.q) - q)))
        @test s.ω === ω
    end
    @info "round trip of arbitrary unit quaternions" max_component_error = worst
    @test worst < 4eps()
end

function spin_sc(m, I, r, v, q, ω)
    bus = mklink(root=true, m=m, dims=(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(r, v; q=q, ang_vel=ω)
    return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=ic, inertia_tensor=SMatrix{3, 3, Float64}(diagm(collect(I))))
end
const QA = normalize(SVector(0.3, -0.5, 0.2, 0.8)); const QB = normalize(SVector(-0.1, 0.4, 0.7, 0.2))
# Spin close to the axis of largest inertia (stable, so discretization errors do not grow exponentially).
const WA = SVector(0.03, 0.02, 0.3); const WB = SVector(0.02, -0.03, 0.25)
spin_scs() = SM.SpacecraftModel[spin_sc(20.0, (2.0, 3.0, 4.0), R1, V1, QA, WA), spin_sc(35.0, (5.0, 6.0, 7.0), R2, V2, QB, WB)]
tcol(df, r, i, name, n) = SVector{n, Float64}(df[r, "sc$(i)_$(name)_$k"] for k in 1:n)

function spin_errors(table, ref, rows)
    dq = 0.0; dw = 0.0
    for r in rows, i in 1:2
        dq = max(dq, quat_angle(tcol(table, r, i, "q", 4), tcol(ref, r, i, "q", 4)))
    end
    return dq
end

@testset "torque-free spinning bodies: ordinary spacecraft vs scene, and negative controls" begin
    effectors = (SM.InverseSquaredGravityModel(),)
    tend = 300.0; every = 25.0
    ref = run_table(engine_config(effectors, tend, every; scs=spin_scs(), orient=true))
    rows = [findmin(abs.(ref.time .- s))[2] for s in every:every:tend]
    function scene_run(dt; states=nothing, integrator=:implicitfast)
        sts = states === nothing ? [body_state_from_initial_condition("chaser", spin_scs()[1].initial_condition),
            body_state_from_initial_condition("target", spin_scs()[2].initial_condition)] : states
        scene = ProximityScene(; mjcf_xml=TWO_BODIES, dt, planet=EARTH, gravity_effectors=effectors, initial_states=sts, integrator, planet_rotation=planet_rotation)
        return run_table(engine_config(effectors, tend, every; scs=spin_scs(), scene, orient=true); extra_callbacks=tick_every(every))
    end
    dts = [0.04, 0.02, 0.01]
    errs = [spin_errors(scene_run(dt), ref, rows) for dt in dts]
    ratios = [errs[i] / errs[i + 1] for i in 1:length(dts) - 1]
    @info "spinning bodies, scene vs ordinary SpaceAGORA over $(tend) s: max attitude error [rad] by scene dt" dts errs ratios
    # MuJoCo's explicit gyroscopic term makes the scene first order in dt (mj_step1/mj_step2 integrate
    # RK4 as Euler), so halving dt halves the error. Measured (rad) 0.0729, 0.0386, 0.0199 at dt = 0.04, 0.02, 0.01
    # and 0.0101 at 0.005; the tolerance at dt = 0.01 is that error with 1.5x margin.
    @test all(r -> 1.7 < r < 2.3, ratios)
    @test errs[end] < 0.03
    # The initial-condition route and the engine hand-off agree: the scene's states in the table are the
    # SpaceAGORA attitude at t = 0 exactly.
    tb = scene_run(0.01)
    @test tcol(tb, 1, 1, "q", 4) == SM.project_unit_quaternion(QA)
    @test tcol(tb, 1, 2, "q", 4) == SM.project_unit_quaternion(QB)
    # Negative controls: a wrong convention is far outside the tolerance. The scene is started from a conjugated
    # quaternion / a world-frame rate, but the engine overwrites scene states from u0, so build the controls
    # from the standalone runner instead and compare with the engine reference directly.
    function standalone_attitudes(states; dt=0.01)
        sc = ProximityScene(; mjcf_xml=TWO_BODIES, dt, planet=EARTH, gravity_effectors=effectors, initial_states=states, planet_rotation=planet_rotation)
        out = Dict{Int, Vector{SVector{4, Float64}}}()
        per = round(Int, every / dt)
        for k in 1:round(Int, tend / dt)
            scene_step!(sc)
            if k % per == 0
                for i in 1:2
                    push!(get!(out, i, SVector{4, Float64}[]), mujoco_to_sa_quaternion(scene_body_state(sc, i).q))
                end
            end
        end
        return out
    end
    function control_error(states)
        att = standalone_attitudes(states)
        return maximum(quat_angle(att[i][k], tcol(ref, rows[k], i, "q", 4)) for i in 1:2, k in eachindex(rows))
    end
    ics = [spin_scs()[1].initial_condition, spin_scs()[2].initial_condition]
    good = [body_state_from_initial_condition(n, ic) for (n, ic) in zip(("chaser", "target"), ics)]
    e_good = control_error(good)
    conj_q(q) = SVector(-q[1], -q[2], -q[3], q[4])
    wrong_conj = [SceneBodyState(s.body, s.r, s.v; q=sa_to_mujoco_quaternion(conj_q(mujoco_to_sa_quaternion(s.q))), ω=s.ω) for s in good]
    R(q) = SM.QuaternionMath.rot(q)
    wrong_world = [SceneBodyState(s.body, s.r, s.v; q=s.q, ω=R(mujoco_to_sa_quaternion(s.q)) * s.ω) for s in good]   # body rate replaced by its inertial image
    wrong_order = [SceneBodyState(s.body, s.r, s.v; q=SVector(s.q[2], s.q[3], s.q[4], s.q[1]), ω=s.ω) for s in good]   # SpaceAGORA order passed to MuJoCo unmapped
    e_conj, e_world, e_order = control_error(wrong_conj), control_error(wrong_world), control_error(wrong_order)
    @info "negative controls: max attitude error [rad]" correct = e_good conjugated = e_conj world_frame_rate = e_world unmapped_order = e_order
    @test e_good < 0.03
    @test e_conj > 0.3 && e_world > 0.3 && e_order > 0.3
end

# --- mixed runs ---------------------------------------------------------------------------------------

# An ordinary spacecraft in a different orbit, and the scene's two bodies after it.
const R3 = SVector(0.0, 7.2e6, 0.0)
const V3_ = sqrt(MU / 7.2e6) .* SVector(-sin(30.0 * π / 180), 0.0, cos(30.0 * π / 180))
ordinary_sc() = rigid_sc(10.0, (1.0, 1.0, 1.0), R3, V3_)
# Loose tolerances with a small dt_max: the step sequence is then set by dt_max alone and does not depend on
# how many other spacecraft share the state vector (the default RMS error norm averages over all components).
loose_tol() = SM.IntegrationTolerances(reltol_orbit=1e-6, abstol_orbit=1e-3, reltol_atmosphere=1e-6, abstol_atmosphere=1e-3,
    reltol_quaternion=1e-6, abstol_quaternion=1e-6, reltol_mass=1e-6, abstol_mass=1e-6,
    reltol_angular_rate=1e-6, abstol_angular_rate=1e-6, dt_max_orbit=1.0, dt_max_atmosphere=1.0)

@testset "mixed run: the ordinary spacecraft is unaffected by the scene to solver round-off" begin
    # Bit-identity to a run of that spacecraft alone is not attainable in SpaceAGORA, scene or not: the solver
    # picks its first step from the RMS norm of the whole state vector, so adding any spacecraft moves that step by
    # ~1e-3 relative and every later step time by ~1e-14 s (about 1e-8 m at orbital speed). The control run shows
    # the same effect with a second ORDINARY spacecraft; the scene must add nothing beyond it.
    effectors = (SM.InverseSquaredJ2GravityModel(),)
    tend = 600.0; every = 10.0
    sc2() = rigid_sc(20.0, (1.0, 0.8, 0.6), R1, V1)
    alone = run_table(engine_config(effectors, tend, every; scs=SM.SpacecraftModel[ordinary_sc()], tol=loose_tol()))
    two = run_table(engine_config(effectors, tend, every; scs=SM.SpacecraftModel[ordinary_sc(), sc2()], tol=loose_tol()))
    scs3 = SM.SpacecraftModel[ordinary_sc(), sc2(), rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)]
    scene = make_scene(effectors; dt=0.05)
    mixed = run_table(engine_config(effectors, tend, every; scs=scs3, scene, owned=(2 => "chaser", 3 => "target"), tol=loose_tol()))
    maxdiff(a, b, f) = maximum(maximum(abs.(a[!, "sc1_$(f)_$k"] .- b[!, "sc1_$(f)_$k"])) for k in 1:3)
    d_two = (maxdiff(alone, two, "pos"), maxdiff(alone, two, "vel"))
    d_mixed = (maxdiff(alone, mixed, "pos"), maxdiff(alone, mixed, "vel"))
    @info "ordinary spacecraft 1: difference from running alone" with_second_ordinary_spacecraft = d_two with_scene = d_mixed
    @test nrow(mixed) == nrow(alone)
    # Measured: 1.49e-8 m and 1.0e-11 m/s with the scene; 8.4e-9 m and 2.9e-11 m/s with a second ordinary spacecraft.
    @test d_mixed[1] < 5e-8 && d_mixed[2] < 3e-11
    @test d_mixed[1] < 10 * max(d_two[1], 1e-9)
    # The scene's bodies are the shadow entries 2 and 3, absolute states in another orbit.
    row = findfirst(==(50.0), mixed.time)
    ref = scene_samples(effectors, 0.05, tend, 50.0)
    @test isapprox(SVector(mixed[row, "sc2_pos_1"], mixed[row, "sc2_pos_2"], mixed[row, "sc2_pos_3"]), ref[1][1]; atol=1e-4)
    @test norm(SVector(mixed[row, "sc2_pos_1"], mixed[row, "sc2_pos_2"], mixed[row, "sc2_pos_3"]) - SVector(mixed[row, "sc1_pos_1"], mixed[row, "sc1_pos_2"], mixed[row, "sc1_pos_3"])) > 1e5
end

# --- Monte Carlo ----------------------------------------------------------------------------------------

@testset "threaded run_monte_carlo gives each sample its own scene and matches serial bit for bit" begin
    effectors = (SM.InverseSquaredGravityModel(),)
    dt = 0.05; tend = 60.0
    template = make_scene(effectors; dt)             # one template shared by every sample and thread
    function sample(seed)
        v1 = V1 + SVector(0.0, 0.0, 1e-3 * seed)
        scs = SM.SpacecraftModel[rigid_sc(20.0, (1.0, 0.8, 0.6), R1, v1), rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)]
        t = run_table(engine_config(effectors, tend, 10.0; scs, scene=template))
        return [t[!, c] for c in names(t) if c != "time"]
    end
    seeds = 1:6
    serial = SpaceAGORA.run_monte_carlo(sample, seeds; threads=1)
    @test isempty(serial.failed)
    nth = min(4, Threads.nthreads())
    if nth > 1
        threaded = SpaceAGORA.run_monte_carlo(sample, seeds; threads=nth)
        @test isempty(threaded.failed)
        @test [s.value for s in threaded.samples] == [s.value for s in serial.samples]      # bit-identical
        @test scene_time(template) == 0.0 && template.n == 0                                # the template was never stepped
    else
        @test_skip false                                                                    # needs julia --threads=N, N > 1
    end
    # different seeds really give different trajectories
    @test serial.samples[1].value != serial.samples[2].value
end

# The chaser carries a hinged link (a body no spacecraft owns); the hinge has no spring, so the link stays put.
const ARM_XML = replace(TWO_BODIES, "</body>\n    <body name=\"target\">" => "<body name=\"link\" pos=\"0 1 0\"><joint type=\"hinge\" axis=\"0 0 1\"/><geom type=\"box\" size=\"0.1 0.5 0.1\" mass=\"1\" contype=\"0\" conaffinity=\"0\"/></body></body>\n    <body name=\"target\">")

@testset "scene_body_pose_save_field: poses of bodies no spacecraft owns" begin
    effectors = (SM.InverseSquaredGravityModel(),)
    dt = 0.05; tend = 20.0; every = 5.0
    scene = ProximityScene(; mjcf_xml=ARM_XML, dt, planet=EARTH, gravity_effectors=effectors, initial_states=states(), planet_rotation=planet_rotation)
    cfg = engine_config(effectors, tend, every; scene)
    fields = vcat(SM.default_save_fields(cfg), [scene_body_pose_save_field("link")])
    tbl = SpaceAGORA.run_simulation(cfg; return_results=true, save_fields=fields, extra_callbacks=tick_every(every)).table
    @test all(c -> c in names(tbl), ["scene_pose_link_$k" for k in 1:7])
    ref = ProximityScene(; mjcf_xml=ARM_XML, dt, planet=EARTH, gravity_effectors=effectors, initial_states=states(), planet_rotation=planet_rotation)
    ok = true
    for (row, t) in enumerate(tbl.time)
        n = round(Int, t / dt)
        while ref.n < n; scene_step!(ref); end
        s = scene_body_state(ref, "link")
        pose = [s.r..., mujoco_to_sa_quaternion(s.q)...]
        ok &= [tbl[row, "scene_pose_link_$k"] for k in 1:7] == pose
    end
    @test ok
    # the link is 1 m from the chaser along y (its frame offset), at the chaser's attitude (identity here)
    last = nrow(tbl)
    @test norm(SVector(tbl[last, "scene_pose_link_1"], tbl[last, "scene_pose_link_2"], tbl[last, "scene_pose_link_3"]) -
               SVector(tbl[last, "sc1_pos_1"], tbl[last, "sc1_pos_2"], tbl[last, "sc1_pos_3"])) ≈ 1.0 atol = 1e-3
    @test_throws ArgumentError SpaceAGORA.run_simulation(cfg; return_results=true, save_fields=[scene_body_pose_save_field("missing")])
end

# --- refusals -----------------------------------------------------------------------------------------------

# A guidance effector that asks for atmosphere-interface events (a continuous event on every spacecraft).
struct EventGuidance end
SpaceAGORA.SimulationLifecycle.requires_atmosphere_events(::EventGuidance) = true

function hinge_sc()
    bus = mklink(root=true, m=10.0, dims=(1.0, 1.0, 1.0))
    panel = mklink(m=2.0, dims=(0.05, 1.0, 0.5), r=(0.0, 1.1, 0.0))
    joint = SM.Joint(bus, SVector(0.0, 0.5, 0.0), panel, SVector(0.0, -0.6, 0.0); joint_type=:hinge, axis=[0, 0, 1])
    return SM.SpacecraftModel(; joints=[joint], links=[bus, panel], root=bus, initial_condition=SM.CartesianInitialCondition(R1, V1))
end

function attached_sc()
    bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
    panel = mklink(m=2.0, dims=(0.05, 1.0, 0.5))
    node = SM.CompliantTopologyNode(:panel; mass_kg=2.0, inertia_body_kg_m2=panel.inertia, position=SVector(0.0, 0.6, 0.0))
    edge = SM.CompliantTopologyEdge(:hold, 0, 1; parent_point_body=SVector(0.0, 0.0, 0.0), child_point_body=SVector(0.0, -0.6, 0.0),
        k_translation_n_m=2e5, c_translation_n_s_m=1e3, k_rotation_n_m_rad=2e4, c_rotation_n_m_s_rad=1e2)
    att = SM.CompliantAttachment(; model=SM.build_compliant_topology([node], [edge]), link=bus, mount_point=(0.0, 0.5, 0.0))
    return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=SM.CartesianInitialCondition(R1, V1),
        inertia_tensor=bus.inertia, attachments=[att])
end

@testset "refusals" begin
    effectors = (SM.InverseSquaredGravityModel(),)
    scene = make_scene(effectors; dt=0.0125)
    cfg(; kw...) = engine_config(effectors, 10.0, 1.0; scene, kw...)
    # The message of the ArgumentError `f` throws (nothing if it throws something else or nothing at all).
    arg_error(f) = try; f(); nothing; catch e; e isa ArgumentError ? sprint(showerror, e) : nothing; end
    refuses(f, rx) = (m = arg_error(f); m !== nothing && occursin(rx, m))
    refused(c, rx; mode=:dp8) = refuses(() -> SE._validate_external_propagators!(c, mode), rx)
    passes(c; mode=:dp8) = SE._validate_external_propagators!(c, mode) === nothing
    @test passes(cfg())                                         # the baseline configuration passes
    @test SpaceAGORA.run_simulation(cfg(); isolate_state=false) === nothing

    # checkpointing and resume
    @test refuses(() -> SpaceAGORA.run_simulation(cfg(extra=(; simulation_settings=SM.SimulationSettings(checkpoint_enabled=true, checkpoint_interval_s=5.0, results=false, results_directory=mktempdir())))), r"checkpoint")
    @test refuses(() -> SpaceAGORA.run_simulation(cfg(extra=(; simulation_settings=SM.SimulationSettings(resume_from_checkpoint=true, results=false, results_directory=mktempdir())))), r"checkpoint")
    # constellation ensembles
    @test refuses(() -> SpaceAGORA.run_constellation_ensemble(cfg()), r"externally propagated")
    # solver routes other than the first-order single-RHS modes
    for mode in (:split_imex, :multirate, :gravity_backbone_split)
        @test refused(cfg(), r"solver mode"; mode)
        @test refuses(() -> SpaceAGORA.run_simulation(cfg(solver=mode)), r"solver mode")
    end
    for mode in (:tsit5, :auto_stiff, :rodas5p, :dp8)
        @test passes(cfg(); mode)
    end
    # the flat constellation route
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "flat") do
        @test refuses(() -> SpaceAGORA.run_simulation(cfg()), r"flat")
    end
    # articulated joints, compliant attachments and the cloth robot arm on a scene-owned spacecraft
    plain = rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)
    @test refuses(() -> SpaceAGORA.run_simulation(cfg(scs=SM.SpacecraftModel[hinge_sc(), plain], orient=true)), r"articulated")
    @test refused(cfg(scs=SM.SpacecraftModel[hinge_sc(), plain], orient=true), r"articulated joints")
    @test refused(cfg(scs=SM.SpacecraftModel[attached_sc(), plain], orient=true), r"compliant attachments")
    arm_model = SM.default_cloth_arm_model(link_lengths_m=(0.9, 0.8, 0.6), link_radii_m=(0.06, 0.05, 0.04), link_masses_kg=(6.0, 4.0, 2.0), mount_offset_body=(0.5, 0.0, 0.6))
    base_pose = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    target = SM.cloth_fk(arm_model, base_pose, [-0.12, -0.08, 0.06]).end_effector_position
    plan = SM.plan_robot_arm_motion(arm_model, base_pose, [0.08, 0.95, -0.85], target; config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=2.0))
    arm = SM.RobotArmControlEffector(plan=plan, spacecraft_idx=1, controller=SM.init_robot_arm_joint_mpc(plan; dt_s=0.1, horizon=6), control_dt_s=0.1)
    @test refused(cfg(extra=(; control_model=SM.ControlModel(control_effectors=(arm,), control_rates=[0.1])), orient=true), r"robot-arm")
    # a control effector that acts on a scene-owned spacecraft, or on every spacecraft
    thr = SM.BaseThrusterModel(thrust=[1.0, 0.0], direction=[1.0, 1.0], Δv=[0.0, 0.0], start_burn_time=[1.0, 1.0], stop_burn_time=[2.0, 2.0], Isp=[300.0, 300.0])
    @test refused(cfg(extra=(; control_model=SM.ControlModel(control_effectors=(thr,), control_rates=[0.1]))), r"Control effector")
    idle = SM.BaseThrusterModel(thrust=[0.0, 0.0], direction=[1.0, 1.0], Δv=[0.0, 0.0], start_burn_time=[1.0, 1.0], stop_burn_time=[2.0, 2.0], Isp=[300.0, 300.0])
    @test passes(cfg(extra=(; control_model=SM.ControlModel(control_effectors=(idle,), control_rates=[0.1]))))   # no thrust on an owned spacecraft
    # continuous events on every spacecraft: orbit-count termination, atmosphere-interface events
    orbits = SM.MissionConfiguration(mission_type=SM.MissionOrbits, keplerian=true, number_of_orbits=1, mission_time=10.0, orientation_sim=false, num_steps_to_save=1000, data_rate=1.0)
    @test refused(cfg(extra=(; mission_configuration=orbits)), r"Orbit-count")
    @test refuses(() -> SpaceAGORA.run_simulation(cfg(extra=(; mission_configuration=orbits))), r"Orbit-count")
    @test refused(cfg(extra=(; guidance_model=SM.GuidanceModel(guidance_effectors=(EventGuidance(),), guidance_rates=[0.0125]))), r"Atmosphere-interface")
    # rates that are not multiples of the scene step
    @test refused(cfg(extra=(; guidance_model=SM.GuidanceModel(guidance_effectors=(thr,), guidance_rates=[0.03]))), r"integer multiple")
    @test passes(cfg(extra=(; guidance_model=SM.GuidanceModel(guidance_effectors=(thr,), guidance_rates=[0.025]))))

    # ownership and construction
    two_scene = ProximitySceneDynamics(scene, 1 => "chaser", 2 => "target")
    bad_index(i) = SM.SimConfig._with_configuration(cfg(); external_propagators=(ProximitySceneDynamics(scene, i => "chaser", 2 => "target"),))
    @test refused(bad_index(3), r"owns spacecraft")                                                         # no such spacecraft
    @test refused(SM.SimConfig._with_configuration(cfg(); external_propagators=(two_scene, two_scene)), r"more than one")   # owned twice
    @test refused(SM.SimConfig._with_configuration(cfg(); external_propagators=(:not_a_propagator,)), r"AbstractExternalPropagator")
    @test_throws ArgumentError ProximitySceneDynamics(scene, 1 => "chaser")             # the target has no spacecraft
    @test_throws ArgumentError ProximitySceneDynamics(scene, 1 => "chaser", 2 => "nope")
    @test_throws ArgumentError ProximitySceneDynamics(scene, 1 => "chaser", 1 => "target")
    @test_throws ArgumentError ProximitySceneDynamics(scene, 1 => "chaser", 2 => "chaser")
    @test_throws ArgumentError ProximitySceneDynamics(scene)
    # mass consistency, when requested
    heavy = SM.SpacecraftModel[rigid_sc(21.0, (1.0, 0.8, 0.6), R1, V1), rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)]
    with_mass(rtol, scs) = SM.SimConfig._with_configuration(cfg(; scs); external_propagators=(ProximitySceneDynamics(scene, 1 => "chaser", 2 => "target"; mass_rtol=rtol),))
    @test refused(with_mass(1e-3, heavy), r"mass")
    @test passes(with_mass(0.1, heavy))
    @test passes(with_mass(1e-9, SM.SpacecraftModel[rigid_sc(20.0, (1.0, 0.8, 0.6), R1, V1), rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)]))
    # a scene with an articulated link: only free-root bodies can carry a spacecraft
    arm_scene = ProximityScene(; mjcf_xml=ARM_XML, dt=0.05, planet=EARTH, gravity_effectors=effectors, initial_states=states(), planet_rotation=planet_rotation)
    @test "link" in scene_body_names(arm_scene)
    @test_throws ArgumentError ProximitySceneDynamics(arm_scene, 1 => "chaser", 2 => "link")
end
