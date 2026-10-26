# Which spacecraft are active in a run resumed from a checkpoint: a spacecraft
# lost to an impact or a touchdown before the checkpoint stays inactive, as in a
# run that was never interrupted (the checkpoint's activity entry).
module CheckpointActivityTests
using Test
using Logging
using Serialization
using StaticArrays
using SpaceAGORA
using OrdinaryDiffEq: DiscreteCallback
const SM = SpaceAGORA.SimulationModel
const CH = SM.ControlHooks
const SE = SpaceAGORA.SimulationEngine
const ZERO3 = SVector(0.0, 0.0, 0.0)
const DIRECTION = SVector(cosd(30.0), 0.0, sind(30.0))

# Navigation, guidance and control that record, per spacecraft, the times they
# tick and the times the engine releases them in place of a tick (an inactive
# spacecraft, with `skip_inactive`), and the activity flags their run-start
# initializers see.
struct Ticks
    ticked::Vector{Vector{Float64}}
    released::Vector{Vector{Float64}}
    started::Vector{Vector{Bool}}
end
Ticks(n) = Ticks([Float64[] for _ in 1:n], [Float64[] for _ in 1:n], Vector{Bool}[])

struct RecordingNavigation <: AbstractNavigationModel
    ticks::Ticks
end
struct RecordingGuidance <: AbstractGuidanceModel
    ticks::Ticks
end
# The control also declares the second spacecraft's touchdown, when given one.
struct RecordingControl <: SM.AbstractControlEffectorModel
    ticks::Ticks
    touchdown::Any
end
SpaceAGORA.update_navigation!(m::RecordingNavigation, u, p, t::Float64, i::Int) = (push!(m.ticks.ticked[i], t); nothing)
SpaceAGORA.update_guidance!(m::RecordingGuidance, u, p, t::Float64, i::Int) = (push!(m.ticks.ticked[i], t); nothing)
SpaceAGORA.update_control!(m::RecordingControl, u, p, t, i) = (push!(m.ticks.ticked[i], t); nothing)
SpaceAGORA.control_wrench(::RecordingControl, x::StateSample, p, t) = (ZERO3, ZERO3)
CH.touchdown_spec(m::RecordingControl, i::Int) = i == 2 ? m.touchdown : nothing
const Recording = Union{RecordingNavigation, RecordingGuidance, RecordingControl}
started!(m::Recording, p) = (push!(m.ticks.started, copy(p.is_active)); nothing)
SpaceAGORA.initialize_navigation!(m::RecordingNavigation, u, p, t) = started!(m, p)
SpaceAGORA.initialize_guidance!(m::RecordingGuidance, u, p, t) = started!(m, p)
SpaceAGORA.initialize_control!(m::RecordingControl, u, p, t) = started!(m, p)
SpaceAGORA.skip_inactive(::Recording) = true      # the 0.4.0 default, declared so the test holds before it
SpaceAGORA.release_command!(m::Recording, t, i) = (push!(m.ticks.released[i], t); nothing)

# Two spacecraft for 10 s, checkpointed every 5 s. The first stays in a 400 km
# orbit. The second comes down at 1 km/s and is lost at about 3 s: an impact at
# 50 km, or a touchdown on terrain 3000 m high at planetocentric 30 degrees (a
# 5 m clearance), whose contacts `touchdowns` records.
function configuration(directory; touchdown::Bool=false, resume::Bool=false, touchdowns=Float64[])
    planet = SM.Earth()
    radius = planet.Rp_e - 2_000.0
    grid = SM.DEMGrid(Float32[3500 3500; 2500 2500], 20, 40, -5, 5; reference_radius_m=radius)
    terrain = SM.DEMTerrainModel([grid]; reference_radius_m=radius)
    spec = touchdown ? (terrain=terrain, reference_radius_m=radius, height_m=5.0,
        on_touchdown=(t, r_p, v_p, i) -> push!(touchdowns, t)) : nothing
    orbit_r = planet.Rp_e + 400e3
    states = ((SVector(orbit_r, 0.0, 0.0), SVector(0.0, sqrt(planet.μ / orbit_r), 0.0)),
              ((planet.Rp_e + (touchdown ? 4_005.0 : 53_000.0)) * DIRECTION, -1000.0 * DIRECTION))
    spacecraft = map(enumerate(states)) do (id, (pos, vel))
        root = SM.Link(root=true, m=500.0, ref_area=1.0)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
            initial_condition=SM.CartesianInitialCondition(pos, vel), id=id)
    end
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false,
            normalize=false, save_csv=false, checkpoint_enabled=true, checkpoint_interval_s=5.0,
            checkpoint_directory=String(directory), resume_from_checkpoint=resume),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=10.0,
            number_of_orbits=1, keplerian=true, orientation_sim=false, num_steps_to_save=20),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
            topography=false, wind=false,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            ephemerides_model=SM.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.0)),
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[spacecraft...], (SM.InverseSquaredJ2GravityModel(),)),
        navigation_model=SM.NavigationModel(navigation_effectors=(RecordingNavigation(Ticks(2)),), navigation_rates=[1.0]),
        guidance_model=SM.GuidanceModel(guidance_effectors=(RecordingGuidance(Ticks(2)),), guidance_rates=[1.0]),
        control_model=SM.ControlModel(control_effectors=(RecordingControl(Ticks(2), spec),), control_rates=[1.0]),
        initial_time=SM.InitialTime(year=2000, month=1, day=1, hour=12, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-10, dt_max_orbit=0.5),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

ticks(args) = [only(args.navigation_model.navigation_effectors).ticks, only(args.guidance_model.guidance_effectors).ticks,
               only(args.control_model.control_effectors).ticks]

# Thrown at the first step past the checkpoint at 5 s, so the run stops there
# as an interrupted run would; unarmed, the callback leaves a run unchanged.
struct Interrupted <: Exception end
interrupter(armed) = DiscreteCallback((u, t, integrator) -> armed[] && t > 5.0, integrator -> throw(Interrupted());
    save_positions=(false, false))
unwrapped(err) = err isa TaskFailedException ? unwrapped(err.task.exception) : err

# A checkpointed run warns that it returns its last segment's solution.
quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))
run!(args, callback) = quiet(() -> run_simulation(args; isolate_state=false, return_solution=true, extra_callbacks=(callback,)))

bits(sol) = [reinterpret(UInt64, x) for u in sol.u for x in collect(Float64, u)]

# The uninterrupted run, and the run interrupted after its checkpoint at 5 s and
# resumed from it, each with its configuration and returned solution (the last
# segment's, from 5 s to 10 s); the interrupted run's directory holds its checkpoint.
function runs(; touchdown::Bool, touchdowns=(Float64[], Float64[], Float64[]))
    armed = Ref(false)
    callback = interrupter(armed)
    straight = configuration(mktempdir(); touchdown, touchdowns=touchdowns[1])
    straight_sol = run!(straight, callback)
    directory = mktempdir()
    interrupted = configuration(directory; touchdown, touchdowns=touchdowns[2])
    armed[] = true
    stopped = try
        run!(interrupted, callback)
        nothing
    catch err
        unwrapped(err)
    finally
        armed[] = false
    end
    resumed = configuration(directory; touchdown, resume=true, touchdowns=touchdowns[3])
    return (; straight, straight_sol, interrupted, stopped, directory, resumed, callback)
end

@testset "A checkpoint keeps which spacecraft are active" begin
    for touchdown in (false, true)
        @testset "$(touchdown ? "touchdown" : "impact") before the checkpoint" begin
            touchdowns = (Float64[], Float64[], Float64[])
            r = runs(; touchdown, touchdowns)
            @test r.stopped isa Interrupted
            checkpoint = SE._load_checkpoint(r.interrupted)
            @test checkpoint.t == 5.0
            @test checkpoint.runtime_state.activity == (version=1, is_active=[true, false])
            @test !hasproperty(checkpoint.runtime_state, :interlinks)
            resumed_sol = run!(r.resumed, r.callback)
            # The second half, bit for bit.
            @test r.straight_sol.t == resumed_sol.t
            @test bits(resumed_sol) == bits(r.straight_sol)
            @test r.straight_sol.prob.p.is_active == resumed_sol.prob.p.is_active == [true, false]
            # The lost spacecraft does not move after its loss.
            @test resumed_sol.u[end].sc[2] == resumed_sol.u[1].sc[2]
            # Every model starts the resumed run with the lost spacecraft inactive,
            # ticks for the spacecraft in orbit as the uninterrupted run does after
            # 5 s, and is released for the lost one in place of a tick.
            for (straight_ticks, resumed_ticks) in zip(ticks(r.straight), ticks(r.resumed))
                @test straight_ticks.started == [[true, true]]
                @test resumed_ticks.started == [[true, false]]
                @test !isempty(resumed_ticks.ticked[1])
                @test resumed_ticks.ticked[1] == filter(>(5.0), straight_ticks.ticked[1])
                @test isempty(resumed_ticks.ticked[2])
                @test !isempty(resumed_ticks.released[2])
                @test resumed_ticks.released[2] == filter(>(5.0), straight_ticks.released[2])
                @test all(<(3.0), straight_ticks.ticked[2])          # ticked until the loss only
            end
            if touchdown
                # One contact before the checkpoint; none again after it.
                @test length(touchdowns[1]) == length(touchdowns[2]) == 1 && only(touchdowns[1]) < 3.0
                @test isempty(touchdowns[3])
            end
        end
    end

    @testset "a checkpoint without the entry resumes as before it" begin
        # An earlier SpaceAGORA wrote no activity entry: the resumed run starts
        # with every spacecraft active, as it did then.
        for touchdown in (false, true)
            r = runs(; touchdown)
            paths = SE._checkpoint_paths(r.interrupted)
            payload = deserialize(paths.data)
            runtime = Base.structdiff(payload.runtime_state, NamedTuple{(:activity,)})
            payload = isempty(runtime) ? Base.structdiff(payload, NamedTuple{(:runtime_state,)}) :
                merge(payload, (runtime_state=runtime,))
            open(io -> serialize(io, payload), paths.data, "w")
            @test !hasproperty(get(SE._load_checkpoint(r.interrupted), :runtime_state, (;)), :activity)
            run!(r.resumed, r.callback)
            @test all(t -> t.started == [[true, true]], ticks(r.resumed))
        end
    end

    @testset "the entry is written with every checkpoint and checked before it is restored" begin
        args = (dynamics_model=(spacecraft=[(;), (;), (;)],), interlink_model=nothing)
        flags = Bool[true, false, true]
        entry = SE._activity_checkpoint_state(flags)
        flags[2] = true
        @test entry == (version=1, is_active=[true, false, true])
        checkpoint(runtime) = (t=5.0, u=nothing, runtime_state=runtime)
        restored = trues(3)
        @test isnothing(SE._restore_activity_checkpoint_state!(args, checkpoint((activity=entry,)), restored))
        @test restored == [true, false, true]
        # No entry: the flags are left as the run set them.
        for legacy in ((t=5.0, u=nothing), checkpoint(nothing), checkpoint((models=nothing,)))
            unchanged = Bool[true, true, false]
            @test isnothing(SE._restore_activity_checkpoint_state!(args, legacy, unchanged))
            @test unchanged == [true, true, false]
        end
        # An interlink entry's flags must agree with it.
        interlinks = (version=1, configuration=nothing, is_active=[true, false, true])
        agreed = trues(3)
        @test isnothing(SE._restore_activity_checkpoint_state!(args, checkpoint((interlinks=interlinks, activity=entry)), agreed))
        @test agreed == [true, false, true]
        for invalid in ((activity=merge(entry, (version=2,)),),
                        (activity=(version=1,),),
                        (activity=merge(entry, (is_active=[true, false],)),),
                        (activity=merge(entry, (is_active=[1, 0, 1],)),),
                        (interlinks=merge(interlinks, (is_active=[true, true, true],)), activity=entry))
            destination = trues(3)
            @test_throws ArgumentError SE._restore_activity_checkpoint_state!(args, checkpoint(invalid), destination)
            @test all(destination)          # a refused entry restores nothing
        end
        @test_throws ArgumentError SE._restore_activity_checkpoint_state!(args, checkpoint((activity=entry,)), trues(2))
    end
end

end # module CheckpointActivityTests
