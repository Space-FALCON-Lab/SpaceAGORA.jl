@testset "free-body agreement with ordinary SpaceAGORA spacecraft: $name" for (name, effectors) in (
        ("point mass", (SM.InverseSquaredGravityModel(),)),
        ("J2", (SM.InverseSquaredJ2GravityModel(),)))
    tend = round(PERIOD / 50) * 50          # one orbit, on a 50 s grid
    every = 50.0
    ref = samples_from_table(engine_reference(effectors, tend, every), every, tend)
    dts = [0.1, 0.05, 0.025, 0.0125]
    diffs = [max_differences(scene_samples(effectors, dt, tend, every), ref) for dt in dts]
    for (dt, d) in zip(dts, diffs)
        @info "$name: scene vs engine over one orbit ($(round(tend)) s), dt = $dt s" max_position_diff_m = d.dr max_velocity_diff_mps = d.dv max_target_position_diff_m = d.dr_target max_relative_position_diff_m = d.drel max_relative_velocity_diff_mps = d.dvrel relative_position = d.dr / 7.0e6
    end
    # Self-convergence by dt halving: the scene's own error (the step holds each wrench for dt) is first
    # order, so successive halvings should roughly halve the difference to the engine.
    ratios = [diffs[i].dr / diffs[i + 1].dr for i in 1:length(dts) - 1]
    @info "$name: position error ratios for dt halving" ratios
    @test all(r -> 1.7 < r < 2.3, ratios)
    @test diffs[end].dr < TOL_POSITION_M
    @test diffs[end].dv < TOL_VELOCITY_MPS
    @test diffs[end].drel < TOL_RELATIVE_POSITION_M
end

# --- runner properties --------------------------------------------------------------------------------

abs_states(sc) = [scene_body_state(sc, i) for i in 1:2]
run_steps!(sc, n) = (for _ in 1:n; scene_step!(sc); end; sc)
same_state(a, b) = all(x.r == y.r && x.v == y.v && x.q == y.q && x.ω == y.ω for (x, y) in zip(abs_states(a), abs_states(b))) &&
    scene_chief(a) == scene_chief(b) && scene_time(a) == scene_time(b)

@testset "deterministic reset" begin
    sc = make_scene((SM.InverseSquaredJ2GravityModel(),); dt=0.05)
    s0 = abs_states(sc)
    run_steps!(sc, 3000)
    first_run = (abs_states(sc), scene_chief(sc))
    @test first_run[1][1].r != s0[1].r
    scene_reset!(sc)
    @test scene_time(sc) == 0.0
    @test all(x.r ≈ y.r && x.v ≈ y.v for (x, y) in zip(abs_states(sc), s0))
    run_steps!(sc, 3000)
    @test (abs_states(sc), scene_chief(sc)) == first_run          # bit-identical
    # a fresh scene reproduces it too, and the initial absolute states are honored
    sc2 = make_scene((SM.InverseSquaredJ2GravityModel(),); dt=0.05)
    run_steps!(sc2, 3000)
    @test same_state(sc, sc2)
    @test maximum(norm(x.r - SVector(R1)) for x in [s0[1]]) < 1e-6
    @test norm(s0[2].r - R2) < 1e-6 && norm(s0[2].v - V2) < 1e-9
end

@testset "chief starts at the mass-weighted COM" begin
    sc = make_scene((SM.InverseSquaredGravityModel(),); dt=0.05)
    R, V = scene_chief(sc)
    @test R ≈ (20 * R1 + 35 * R2) / 55 atol = 1e-6
    @test V ≈ (20 * V1 + 35 * V2) / 55 atol = 1e-9
end

@testset "state get/set round trip" begin
    sc = make_scene((SM.InverseSquaredJ2GravityModel(),); dt=0.05)
    run_steps!(sc, 700)
    st = scene_state(sc)
    @test st.n == 700
    run_steps!(sc, 900)
    ref = (abs_states(sc), scene_chief(sc), scene_time(sc))
    scene_set_state!(sc, st)
    @test scene_time(sc) == 700 * 0.05
    run_steps!(sc, 900)
    @test (abs_states(sc), scene_chief(sc), scene_time(sc)) == ref
    # copy() gives an independent scene (own mjModel and mjData) at the same state
    scene_set_state!(sc, st)
    cp = copy(sc)
    run_steps!(cp, 900)
    @test sc.n == 700 && cp.n == 1600
    @test (abs_states(cp), scene_chief(cp), scene_time(cp)) == ref
    @test_throws DimensionMismatch scene_set_state!(sc, SceneState(0, st.R, st.V, zeros(3)))
end

@testset "configuration is validated" begin
    kw = (; mjcf_xml=TWO_BODIES, planet=EARTH, gravity_effectors=(SM.InverseSquaredGravityModel(),), initial_states=states())
    @test_throws UndefKeywordError ProximityScene(; kw...)                              # dt is required
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.0)
    @test_throws ArgumentError ProximityScene(; kw..., dt=NaN)
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.05, integrator=:midpoint)
    # RK4 would silently run as Euler on the mj_step1/mj_step2 path, so it is refused
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.05, integrator=:rk4)
    @test_throws ArgumentError ProximityScene(; kw..., mjcf_path="x.xml", dt=0.05)    # both sources
    @test_throws ArgumentError ProximityScene(; planet=EARTH, gravity_effectors=kw.gravity_effectors, initial_states=states(), dt=0.05)
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.05, initial_states=states()[1:1])   # missing body state
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.05, initial_states=[states()[1], states()[1]])
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.05, gravity_effectors=(SM.InverseSquaredJ2GravityModel(),))  # no planet_rotation
    @test_throws ArgumentError ProximityScene(; kw..., dt=0.05, gravity_effectors=())
    # the default integrator is implicitfast
    @test make_scene((SM.InverseSquaredGravityModel(),); dt=0.05).integrator === :implicitfast
    sc = make_scene((SM.InverseSquaredGravityModel(),); dt=0.05)
    @test Binding.timestep(sc.model) == 0.05 && Binding.gravity(sc.model) == (0.0, 0.0, 0.0)
    @test Binding.integrator(sc.model) == Int(Binding.INT_IMPLICITFAST)
    # SpaceAGORA initial conditions map to scene states, attitude and rate included (the convention is
    # derived and tested in engine_tests.jl)
    ic_ok = SM.CartesianInitialCondition(R1, V1)
    @test body_state_from_initial_condition("chaser", ic_ok).r == R1
    @test body_state_from_initial_condition("chaser", ic_ok).q == SVector(1.0, 0.0, 0.0, 0.0)
    ic_rot = SM.CartesianInitialCondition(R1, V1; q=SVector(0.0, 0.0, 0.6, 0.8), ang_vel=SVector(0.0, 0.0, 0.1))
    st = body_state_from_initial_condition("chaser", ic_rot)
    @test st.q ≈ SVector(0.8, 0.0, 0.0, 0.6) && st.ω == SVector(0.0, 0.0, 0.1)
end

@testset "feedback acceleration keeps the chief on a thrusted body" begin
    xml = """<mujoco><worldbody><body name="solo"><freejoint/><inertial pos="0 0 0" mass="20" diaginertia="1 1 1"/></body></worldbody></mujoco>"""
    mk() = ProximityScene(; mjcf_xml=xml, dt=0.05, planet=EARTH, gravity_effectors=(SM.InverseSquaredGravityModel(),),
        initial_states=[SceneBodyState("solo", R1, V1)])
    coast, burn = mk(), mk()
    F = reshape([0.5, 0.0, 0.0], 3, 1)
    for _ in 1:2000
        scene_step!(coast); scene_step!(burn; external_forces=F)
    end
    # the single body is the chief: no relative offset develops, thrust acceleration F/m moves both
    @test abs(burn.xipos[1, 2]) + abs(burn.xipos[2, 2]) + abs(burn.xipos[3, 2]) < 1e-9
    dr = scene_body_state(burn, 1).r - scene_body_state(coast, 1).r
    @test dr[1] ≈ 0.5 * (0.5 / 20) * 100.0^2 rtol = 1e-2
    @test scene_chief(burn)[1] ≈ scene_body_state(burn, 1).r
end

@testset "native memory: repeated make/free neither leaks nor crashes" begin
    rss_mb() = parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 2^20
    cycle() = for _ in 1:150
        sc = make_scene((SM.InverseSquaredGravityModel(),); dt=0.05)
        run_steps!(sc, 3)
        cp = copy(sc)
        run_steps!(cp, 3)
    end
    cycle(); GC.gc(true); GC.gc(true)
    before = rss_mb()
    for _ in 1:4
        cycle(); GC.gc(true); GC.gc(true)
    end
    growth = rss_mb() - before
    @info "RSS growth over 1200 scene make/free cycles (MB)" growth
    @test growth < 60
end

@testset "MuJoCo instability is an error, not a silent reset" begin
    sc = make_scene((SM.InverseSquaredJ2GravityModel(),); dt=0.05)
    run_steps!(sc, 3)
    forces = zeros(3, 2); forces[1, 1] = Inf       # drives qacc to Inf; MuJoCo would reset the data and go on
    err = try scene_step!(sc; external_forces=forces); nothing catch e; e end
    @test err isa ErrorException
    @test occursin("unstable", err.msg) && occursin("mjWARN_BADQ", err.msg) && occursin("step 4", err.msg)
    # a huge but finite force is caught the same way
    sc = make_scene((SM.InverseSquaredJ2GravityModel(),); dt=0.05)
    forces = fill(1e300, 3, 2)
    @test_throws ErrorException scene_step!(sc; external_forces=forces)
end
