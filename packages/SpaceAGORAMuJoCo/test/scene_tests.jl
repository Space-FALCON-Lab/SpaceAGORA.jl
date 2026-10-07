using LinearAlgebra
using SpaceAGORAMuJoCo: Binding
using StaticArrays
using DataFrames
import SpaceAGORA
using SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const MU = 3.986004418e14

# Two free rigid boxes, body frames at their centers of mass, no contact, no joints.
const TWO_BODIES = """
<mujoco>
  <compiler angle="degree"/>
  <option timestep="0.05"/>
  <worldbody>
    <body name="chaser"><freejoint/>
      <inertial pos="0 0 0" mass="20" diaginertia="2 3 4"/>
      <geom type="box" size="0.5 0.4 0.3" contype="0" conaffinity="0" mass="0"/>
    </body>
    <body name="target"><freejoint/>
      <inertial pos="0 0 0" mass="35" diaginertia="5 6 7"/>
      <geom type="box" size="0.6 0.6 0.6" contype="0" conaffinity="0" mass="0"/>
    </body>
  </worldbody>
</mujoco>"""

const EARTH = SM.make_no_gram_planet(:earth)
const EPH = SM.SimpleEphemeridesModel()
const T0 = SM.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)
const ET0 = SM.ephemerides_time_seconds(T0, EPH)
planet_rotation(t) = SM.planet_frame_lpi(EARTH, ET0 + t, EPH)

# chaser on a 51.6 deg circular orbit at r = 7000 km; the target is ~50 m away with a small velocity offset.
const R1 = SVector(7.0e6, 0.0, 0.0)
const V1 = sqrt(MU / 7.0e6) .* SVector(0.0, cosd(51.6), sind(51.6))
const R2 = R1 + SVector(30.0, -40.0, 15.0)
const V2 = V1 + SVector(0.05, 0.02, -0.03)
# Tolerances at dt = 0.0125 s, set from the measured first-order convergence (errors scale as dt; measured
# 1.70e-3 m, 2.17e-6 m/s and 2.68e-3 m over one orbit, point mass and J2), with 1.5x margin.
const TOL_POSITION_M = 2.5e-3
const TOL_VELOCITY_MPS = 3.3e-6
const TOL_RELATIVE_POSITION_M = 4.0e-3
const PERIOD = 2π * sqrt(7.0e6^3 / MU)

states() = [SceneBodyState("chaser", R1, V1), SceneBodyState("target", R2, V2)]
make_scene(effectors; dt, kwargs...) = ProximityScene(; mjcf_xml=TWO_BODIES, dt, planet=EARTH, gravity_effectors=effectors,
    initial_states=states(), planet_rotation=planet_rotation, kwargs...)

# --- the same two bodies as ordinary SpaceAGORA spacecraft ---------------------------------------------

mklink(; root=false, m, dims) = SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(0, 0, 0), q=MVector{4, Float64}(0, 0, 0, 1))
function rigid_sc(m, dims, r, v)
    bus = mklink(root=true, m=m, dims=dims)
    ic = SM.CartesianInitialCondition(r, v)
    return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=ic, inertia_tensor=bus.inertia)
end
function engine_reference(effectors, tend, data_rate)
    tol = SM.IntegrationTolerances(reltol_orbit=1e-12, abstol_orbit=1e-9, reltol_atmosphere=1e-12, abstol_atmosphere=1e-9,
        reltol_quaternion=1e-12, abstol_quaternion=1e-9, reltol_mass=1e-12, abstol_mass=1e-9,
        reltol_angular_rate=1e-12, abstol_angular_rate=1e-9, dt_max_orbit=0.5, dt_max_atmosphere=0.5)
    scs = SM.SpacecraftModel[rigid_sc(20.0, (1.0, 0.8, 0.6), R1, V1), rigid_sc(35.0, (1.2, 1.2, 1.2), R2, V2)]
    base = make_example_config(planet=EARTH, spacecraft=scs[1], mission_time=tend, initial_time=T0,
        dynamic_effectors=effectors, density_model=SM.NoAtmosphereModel(), ephemerides_model=EPH,
        orientation_sim=false, keplerian=true, verbose=false, results=false, results_directory=mktempdir(),
        solver_config=SM.SolverConfig(solver_mode=:dp8))
    args = SM.SimConfig._with_configuration(base;
        integration_tolerances=tol,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=tend, orientation_sim=false, num_steps_to_save=100000, data_rate=data_rate),
        dynamics_model=SM.DynamicsModel(scs, effectors))
    return SpaceAGORA.run_simulation(args; return_results=true).table
end

# Run the scene to `tend`, sampling absolute states every `every` seconds (a whole number of steps).
function scene_samples(effectors, dt, tend, every)
    sc = make_scene(effectors; dt)
    per = round(Int, every / dt); @assert per * dt ≈ every
    nsteps = round(Int, tend / dt)
    out = Vector{NTuple{4, SVector{3, Float64}}}()
    for k in 1:nsteps
        scene_step!(sc)
        if k % per == 0
            a = scene_body_state(sc, "chaser"); b = scene_body_state(sc, "target")
            push!(out, (a.r, a.v, b.r, b.v))
        end
    end
    return out
end

col(df, k) = [df[!, "sc1_pos_$k"], df[!, "sc1_vel_$k"], df[!, "sc2_pos_$k"], df[!, "sc2_vel_$k"]]
function samples_from_table(df, every, tend)
    t = Float64.(df[!, :time])
    idx = [findmin(abs.(t .- s))[2] for s in every:every:tend]
    all(isapprox.(t[idx], collect(every:every:tend); atol=1e-6)) || error("engine output times do not line up")
    c = [col(df, k) for k in 1:3]
    return [(SVector(c[1][1][i], c[2][1][i], c[3][1][i]), SVector(c[1][2][i], c[2][2][i], c[3][2][i]),
             SVector(c[1][3][i], c[2][3][i], c[3][3][i]), SVector(c[1][4][i], c[2][4][i], c[3][4][i])) for i in idx]
end

function max_differences(a, b)
    dr = maximum(norm(x[1] - y[1]) for (x, y) in zip(a, b)); dv = maximum(norm(x[2] - y[2]) for (x, y) in zip(a, b))
    dr2 = maximum(norm(x[3] - y[3]) for (x, y) in zip(a, b))
    drel = maximum(norm((x[3] - x[1]) - (y[3] - y[1])) for (x, y) in zip(a, b))
    dvrel = maximum(norm((x[4] - x[2]) - (y[4] - y[2])) for (x, y) in zip(a, b))
    return (dr = dr, dv = dv, dr_target = dr2, drel = drel, dvrel = dvrel)
end

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
    # attitude/rate mapping from SpaceAGORA initial conditions is refused instead of dropped
    ic_ok = SM.CartesianInitialCondition(R1, V1)
    @test body_state_from_initial_condition("chaser", ic_ok).r == R1
    ic_rot = SM.CartesianInitialCondition(R1, V1; ang_vel=SVector(0.0, 0.0, 0.1))
    @test_throws ArgumentError body_state_from_initial_condition("chaser", ic_rot)
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
