module ExternalPropagationTests

# The generic externally-propagated-spacecraft versioned service contract with a trivial,
# MuJoCo-free owner: a spacecraft on a circular orbit about the planet with a constant spin about z, advanced in
# whole steps of `dt` and followed by a shadow entry in the engine state.

using Test
using LinearAlgebra
using StaticArrays
using DataFrames
using SpaceAGORA
using SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const EP = SpaceAGORA.ExternalPropagationServices
const SE = SpaceAGORA.SimulationEngine

const MU = 3.986004418e14
const RADIUS = 7.0e6
const OMEGA_ORBIT = sqrt(MU / RADIUS^3)
const SPIN = 0.05                                    # rad/s about the body z axis
const T0 = SM.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)

# Analytic owner state at time t: circular orbit in the x-y plane, spin about z (body to inertial q, scalar-last).
orbit_pos(t) = RADIUS * SVector(cos(OMEGA_ORBIT * t), sin(OMEGA_ORBIT * t), 0.0)
orbit_vel(t) = RADIUS * OMEGA_ORBIT * SVector(-sin(OMEGA_ORBIT * t), cos(OMEGA_ORBIT * t), 0.0)
spin_quat(t) = SVector(0.0, 0.0, sin(SPIN * t / 2), cos(SPIN * t / 2))

struct CircleOwner <: EP.AbstractExternalPropagator
    owned::Vector{Int}
    dt::Float64
end
mutable struct CircleRuntime
    dt::Float64
    n::Int
    sync_times::Vector{Float64}
    state_calls::Vector{Tuple{Float64, Float64}}     # (engine time, owner time) at each external_state call
    prepared_from::Any
end
EP.external_spacecraft(o::CircleOwner) = copy(o.owned)
EP.external_step(o::CircleOwner) = o.dt
function EP.external_prepare(o::CircleOwner, args, u0)
    return CircleRuntime(o.dt, 0, Float64[], Tuple{Float64, Float64}[], EP.external_initial_state(u0, o.owned[1]).pos)
end
function EP.external_sync!(rt::CircleRuntime, t::Float64)
    push!(rt.sync_times, t)
    target = floor(Int, t / rt.dt + 1e-6)
    steps = target - rt.n
    rt.n = max(rt.n, target)
    return max(steps, 0)
end
EP.external_time(rt::CircleRuntime) = rt.n * rt.dt
function EP.external_state(rt::CircleRuntime, k::Int, t::Float64)
    push!(rt.state_calls, (t, rt.n * rt.dt))
    return (pos = orbit_pos(t), vel = orbit_vel(t), q = spin_quat(t), ω = SVector(0.0, 0.0, SPIN))
end
EP.external_acceleration(rt::CircleRuntime, k::Int, t::Float64) = -OMEGA_ORBIT^2 * orbit_pos(t)

mklink(; root=false, m, dims) = SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(0, 0, 0), q=MVector{4, Float64}(0, 0, 0, 1))
function sc_at(pos, vel; m=10.0, q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, 0.0))
    bus = mklink(root=true, m=m, dims=(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(pos, vel; q=q, ang_vel=ang_vel)
    return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus, initial_condition=ic, inertia_tensor=bus.inertia)
end

# Spacecraft 1 is ordinary on another circular orbit; spacecraft 2 sits on the owner's circle.
function config(owner; tend=300.0, data_rate=50.0, orient=true, solver=:dp8, extra=(;))
    ordinary = sc_at(SVector(0.0, 7.4e6, 0.0), sqrt(MU / 7.4e6) * SVector(-1.0, 0.0, 0.0))
    owned = sc_at(orbit_pos(0.0), orbit_vel(0.0); q=spin_quat(0.0), ang_vel=SVector(0.0, 0.0, SPIN))
    scs = SM.SpacecraftModel[ordinary, owned]
    effectors = (SM.InverseSquaredGravityModel(),)
    base = make_example_config(planet=SM.make_no_gram_planet(:earth), spacecraft=ordinary, mission_time=tend, initial_time=T0,
        dynamic_effectors=effectors, density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=orient, keplerian=true, verbose=false, results=false, results_directory=mktempdir(),
        solver_config=SM.SolverConfig(solver_mode=solver))
    updated = SM.SimConfig._with_configuration(base;
        mission_configuration=SM.MissionConfiguration(mission_type=base.mission_configuration.mission_type, keplerian=true,
            number_of_orbits=1, mission_time=tend, orientation_sim=orient, num_steps_to_save=100000, data_rate=data_rate),
        dynamics_model=SM.DynamicsModel(scs, effectors), extra...)
    owners = get(extra, :external_propagators, owner === nothing ? () : (owner,))
    return EP.with_external_propagators(updated, owners)
end

arg_error(f) = try; f(); nothing; catch e; e isa ArgumentError ? sprint(showerror, e) : nothing; end
refuses(f, rx) = (m = arg_error(f); m !== nothing && occursin(rx, m))

@testset "versioned external propagation service boundary" begin
    @test EP.CONTRACT_VERSION == v"1.0.0"
    @test EP.require_version(v"1.0.0") === nothing
    @test_throws EP.CompatibilityError EP.require_version(v"2.0.0")
    mismatch = try; EP.require_version(v"2.0.0"); catch error; sprint(showerror, error); end
    @test occursin("2.0.0", mismatch) && occursin("1.0.0", mismatch)
    @test EP.AbstractExternalPropagator === SM.ExternalPropagation.AbstractExternalPropagator
    for name in EP.IMPLEMENTED_FUNCTIONS
        @test getfield(EP, name) === getfield(SM.ExternalPropagation, name)
    end
    @test EP.sample_planet_frame_with_lpi === SE.sample_planet_frame_with_lpi
    @test EP.gravity_backbone_relative_acceleration_ii === SM.gravity_backbone_relative_acceleration_ii
    @test isempty(intersect(Set(EP.IMPLEMENTED_FUNCTIONS), Set(EP.CONSUMED_SERVICES)))
    metadata = Set((:ExternalPropagationServices, :CONTRACT_VERSION, :CompatibilityError,
        :require_version, :IMPLEMENTED_FUNCTIONS, :CONSUMED_SERVICES))
    @test Set(names(EP)) == union(metadata, Set(EP.IMPLEMENTED_FUNCTIONS), Set(EP.CONSUMED_SERVICES))
    for name in union(metadata, Set(EP.IMPLEMENTED_FUNCTIONS), Set(EP.CONSUMED_SERVICES))
        @test Base.Docs.doc(Base.Docs.Binding(EP, name)) !== nothing
    end

    original = config(nothing)
    owners = (CircleOwner([2], 0.5),)
    updated = EP.with_external_propagators(original, owners)
    @test original.external_propagators === ()
    @test updated.external_propagators === owners
    for field in fieldnames(typeof(original))
        field === :external_propagators && continue
        @test getfield(updated, field) === getfield(original, field)
    end
    @test EP.with_external_propagators(updated, ()).external_propagators === ()
    @test EP.validate_external_propagators(updated, :dp8) === nothing

    # Reading initial state does not expose mutable engine storage, and disabled
    # attitude has the same identity/zero defaults as an ordinary spacecraft.
    translation = (pos=[1.0, 2.0, 3.0], vel=[4.0, 5.0, 6.0])
    rotation = (pos=[7.0, 8.0, 9.0], vel=[10.0, 11.0, 12.0],
        q=[0.0, 0.0, sin(0.2), cos(0.2)], ω=[0.1, 0.2, 0.3])
    initial = (sc=[translation, rotation],)
    state1 = EP.external_initial_state(initial, 1)
    state2 = EP.external_initial_state(initial, 2)
    @test state1 == (pos=SVector(1.0, 2.0, 3.0), vel=SVector(4.0, 5.0, 6.0),
        q=SVector(0.0, 0.0, 0.0, 1.0), ω=zero(SVector{3, Float64}))
    @test state2.pos == rotation.pos && state2.vel == rotation.vel
    @test state2.q == rotation.q && state2.ω == rotation.ω
    @test all(value -> value isa SVector, values(state2))
    rotation.pos[1] = 100.0
    rotation.q[1] = 0.5
    @test state2.pos[1] == 7.0 && state2.q[1] == 0.0
    @test_throws BoundsError EP.external_initial_state(initial, 3)

    # Consumers receive the run's objects in owner order, not template copies.
    runtimes = Any[CircleRuntime(0.5, 0, Float64[], Tuple{Float64, Float64}[], nothing),
        CircleRuntime(1.0, 0, Float64[], Tuple{Float64, Float64}[], nothing)]
    integrator = (p=(shared_buffers=(external_runtimes=runtimes,),),)
    @test collect(EP.external_runtimes(integrator)) == runtimes
    @test first(EP.external_runtimes(integrator)) === first(runtimes)
    @test isempty(EP.external_runtimes((p=(shared_buffers=(external_runtimes=Any[],),),)))
end

@testset "public planet frame sampling preserves rotating-frame velocity" begin
    planet = make_no_gram_planet(:earth)
    l_pi = @SMatrix [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0]
    pos_pp = SVector(planet.Rp_e + 500.0, 0.0, 0.0)
    vel_pp = SVector(12.0, 34.0, 56.0)
    state = StateSample(l_pi' * pos_pp, l_pi' * (vel_pp + cross(planet.ω, pos_pp)), 1.0)
    frame = EP.sample_planet_frame_with_lpi(state, planet, l_pi)
    @test frame.l_pi === l_pi
    @test frame.pos_pp == pos_pp
    @test frame.vel_pp ≈ vel_pp atol=1e-12 rtol=0
    @test frame.alt_m ≈ 500.0 atol=1e-6 rtol=0
    @test frame.lat_rad ≈ 0.0 atol=1e-12
    @test frame.lon_rad ≈ 0.0 atol=1e-12
end

@testset "hook defaults throw Not implemented" begin
    struct Bare <: EP.AbstractExternalPropagator end
    bare = Bare()
    for f in (EP.external_spacecraft, EP.external_step)
        @test_throws ErrorException f(bare)
    end
    @test_throws ErrorException EP.external_prepare(bare, nothing, nothing)
    @test_throws ErrorException EP.external_sync!(nothing, 1.0)
    @test_throws ErrorException EP.external_time(nothing)
    @test_throws ErrorException EP.external_state(nothing, 1, 1.0)
    @test_throws ErrorException EP.external_acceleration(nothing, 1, 1.0)
    @test EP.external_preflight(bare, nothing) === nothing        # the one hook with a usable default
    msg = try; EP.external_step(bare); catch e; sprint(showerror, e); end
    @test occursin("Not implemented", msg)
    @test config(nothing).external_propagators === ()
end

@testset "shadow entry follows the owner; the ordinary spacecraft is untouched" begin
    dt = 0.5
    tend = 300.0
    owner = CircleOwner([2], dt)
    cfg = config(owner; tend)
    res = SpaceAGORA.run_simulation(cfg; return_results=true, isolate_state=false)
    tbl = res.table
    # shadow position and attitude follow the analytic circle and spin
    perr = maximum(norm(SVector(tbl[i, "sc2_pos_1"], tbl[i, "sc2_pos_2"], tbl[i, "sc2_pos_3"]) - orbit_pos(tbl.time[i])) for i in 1:nrow(tbl))
    qerr = maximum(norm(SVector(tbl[i, "sc2_q_1"], tbl[i, "sc2_q_2"], tbl[i, "sc2_q_3"], tbl[i, "sc2_q_4"]) - spin_quat(tbl.time[i])) for i in 1:nrow(tbl))
    @info "shadow vs analytic owner" position_error_m = perr quaternion_error = qerr
    @test perr < 1e-6                                   # measured 1.0e-9 m
    @test qerr < 2e-8                                   # measured 3.8e-9: the quaternion integration tolerance
    # the ordinary spacecraft (index 1) matches a run in which nothing is owned, to solver round-off
    alone = SpaceAGORA.run_simulation(config(nothing; tend); return_results=true).table
    d = maximum(abs.(tbl[!, "sc1_pos_$k"] .- alone[!, "sc1_pos_$k"]) |> maximum for k in 1:3)
    @info "ordinary spacecraft: owned run vs unowned run" max_position_difference_m = d
    @test d < 1e-6
    # hook protocol: with isolate_state=false the run's runtime is the one built from the initial state of the
    # owned spacecraft; it is synced after accepted steps with nondecreasing times, and every state request
    # is for a time within one owner step after the owner's own time
    rt = nothing
    probe = SM.SimulationCallbacks.PeriodicCallback(integrator -> (rt = first(EP.external_runtimes(integrator))), 100.0)
    SpaceAGORA.run_simulation(config(owner; tend); extra_callbacks=(probe,), isolate_state=false)
    @test rt isa CircleRuntime
    @test rt.prepared_from == orbit_pos(0.0)
    @test issorted(rt.sync_times) && !isempty(rt.sync_times)
    @test rt.n * dt == floor(tend / dt + 1e-6) * dt
    @test all(0 <= te - to < dt + 1e-9 for (te, to) in rt.state_calls)
end

@testset "refusals of the generic hook" begin
    owner = CircleOwner([2], 0.5)
    @test EP.validate_external_propagators(config(owner), :dp8) === nothing
    # checkpoint and resume
    for settings in (SM.SimulationSettings(checkpoint_enabled=true, checkpoint_interval_s=50.0, results=false, results_directory=mktempdir()),
            SM.SimulationSettings(resume_from_checkpoint=true, results=false, results_directory=mktempdir()))
        @test refuses(() -> SpaceAGORA.run_simulation(config(owner; extra=(; simulation_settings=settings))), r"checkpoint")
    end
    # constellation ensemble
    @test refuses(() -> SpaceAGORA.run_constellation_ensemble(config(owner)), r"externally propagated")
    # solver routes
    for mode in (:split_imex, :multirate, :gravity_backbone_split)
        @test refuses(() -> EP.validate_external_propagators(config(owner), mode), r"solver mode")
    end
    for mode in (:tsit5, :auto_stiff, :rodas5p, :dp8)
        @test EP.validate_external_propagators(config(owner), mode) === nothing
    end
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "flat") do
        @test refuses(() -> SpaceAGORA.run_simulation(config(owner)), r"flat")
    end
    # ownership
    @test refuses(() -> EP.validate_external_propagators(config(CircleOwner([3], 0.5)), :dp8), r"owns spacecraft")
    @test refuses(() -> EP.validate_external_propagators(config(CircleOwner([2, 2], 0.5)), :dp8), r"more than one|twice")
    @test refuses(() -> EP.validate_external_propagators(config(CircleOwner(Int[], 0.5)), :dp8), r"owns no spacecraft")
    @test refuses(() -> EP.validate_external_propagators(config(CircleOwner([2], -1.0)), :dp8), r"step")
    @test refuses(() -> EP.validate_external_propagators(config(owner; extra=(; external_propagators=(1,))), :dp8), r"AbstractExternalPropagator")
    # a control effector that acts on an owned (or every) spacecraft; one with no thrust on it is fine
    thr(t1) = SM.BaseThrusterModel(thrust=[0.0, t1], direction=[1.0, 1.0], Δv=[0.0, 0.0], start_burn_time=[1.0, 1.0], stop_burn_time=[2.0, 2.0], Isp=[300.0, 300.0])
    ctl(m) = (; control_model=SM.ControlModel(control_effectors=(m,), control_rates=[0.5]))
    @test refuses(() -> EP.validate_external_propagators(config(owner; extra=ctl(thr(1.0))), :dp8), r"Control effector")
    @test EP.validate_external_propagators(config(owner; extra=ctl(thr(0.0))), :dp8) === nothing
    # guidance and navigation rates must be multiples of the owner step
    gdn(rate) = (; guidance_model=SM.GuidanceModel(guidance_effectors=(thr(0.0),), guidance_rates=[rate]))
    @test refuses(() -> EP.validate_external_propagators(config(owner; extra=gdn(0.7)), :dp8), r"integer multiple")
    @test EP.validate_external_propagators(config(owner; extra=gdn(1.5)), :dp8) === nothing
    # orbit-count termination is a continuous event on every spacecraft
    orbits = SM.MissionConfiguration(mission_type=SM.MissionOrbits, keplerian=true, number_of_orbits=1, mission_time=100.0, orientation_sim=true, num_steps_to_save=100, data_rate=10.0)
    @test refuses(() -> EP.validate_external_propagators(config(owner; extra=(; mission_configuration=orbits)), :dp8), r"Orbit-count")
end

end # module ExternalPropagationTests
