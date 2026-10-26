# The maneuver guidance models declared on spacecraft (SpaceAGORA 0.7.0): each
# copy is bound to its spacecraft and acts on it only, beside one burn scheduler for
# every spacecraft. Copies with equal settings give the bits of one model at
# configuration level; copies with their own settings each follow their own. A
# bound copy that still acts on every spacecraft warns, and an ensemble refuses a
# configuration-level model bound to run positions.
module BoundManeuverGuidanceTests
using Test, Logging, StaticArrays, LinearAlgebra, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance
const L = SpaceAGORA.SimulationLifecycle
const Z3 = SVector(0.0, 0.0, 0.0)
const PLANET = SM.make_no_gram_planet(:earth)

raise(; h_p=500e3) = SM.ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=PLANET.Rp_e + 3000e3,
    target_periapsis_altitude_m=h_p)
campaign(orbits, Δv) = SM.AerobrakingCampaignPropulsiveManeuverGuidanceModel(maneuver_orbit_number=orbits,
    maneuver_Δv=Δv)

# Three spacecraft on 400 km by 2000 km orbits with one burn scheduler for all of
# them, and guidance at configuration level (`shared`) or declared on each
# spacecraft (`declared[k]` on spacecraft k), at a 30 s rate. The scheduler plans
# at its control ticks after due guidance: a periapsis raise plans from inside
# its window before apoapsis (`:raise` starts
# the spacecraft there), a campaign maneuver after exits of an interface above
# the periapsis (`:campaign`; no atmosphere, so no drag).
function maneuver_args(scenario::Symbol=:raise; shared=(), declared=((), (), ()), mission_time=14000.0, settings=(;))
    sats = map(1:3) do k
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=PLANET.Rp_e + 2000e3, rp=PLANET.Rp_e + 400e3, i=30.0, ω=0.0, Ω=40.0k,
            ν=scenario === :raise ? 150.0 + 5.0k : 10.0 + 30.0k)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=50.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=k,
            guidance=SM.GuidanceModel(Tuple(declared[k]), fill(30.0, length(declared[k]))))
    end
    scheduler = SM.ThrusterModels.BaseThrusterModel(thrust=fill(50.0, 3), direction=zeros(3), Δv=zeros(3),
        start_burn_time=fill(-1.0, 3), stop_burn_time=fill(-1.0, 3), Isp=fill(300.0, 3))
    env = SM.EnvironmentModel(planet=PLANET, EI=scenario === :raise ? 120.0 : 500.0,
        density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(; results=false, verbose=false, generate_plots=false, normalize=false,
            settings...),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=1000),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(Tuple(shared), fill(30.0, length(shared))),
        navigation_model=SM.NavigationModel((), Float64[]),
        control_model=SM.ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

solve(args) = run_simulation(args; return_solution=true)

# The saved states, one per saved time (the last at a time saved more than once:
# each periodic callback that fires at a time saves there).
function by_time(sol)
    keep = [i for i in eachindex(sol.t) if i == lastindex(sol.t) || sol.t[i] != sol.t[i + 1]]
    return sol.t[keep], [collect(sol.u[i]) for i in keep]
end
same_bits(a, b) = by_time(a)[1] == by_time(b)[1] && all(isequal.(by_time(a)[2], by_time(b)[2]))

# Spacecraft k's osculating periapsis and semimajor axis at the end of the run.
function final_orbit(sol, k)
    pos, vel = SVector{3}(sol.u[end].sc[k].pos), SVector{3}(sol.u[end].sc[k].vel)
    a = inv(2 / norm(pos) - dot(vel, vel) / PLANET.μ)
    e = norm(cross(vel, cross(pos, vel)) / PLANET.μ - pos / norm(pos))
    return (periapsis=a * (1 - e), a=a)
end

# The log records `f()` emits at warning level that mention `needle`.
function warnings(f, needle)
    logger = Test.TestLogger(min_level=Logging.Warn)
    Logging.with_logger(f, logger)
    return [r for r in logger.logs if r.level == Logging.Warn && occursin(needle, string(r.message))]
end
const BOUND_ALL = "and bound to it by bind_spacecraft, but its spacecraft_scope is AllSpacecraft()"

# A user model whose bound copy keeps the default scope, and one whose copy does not.
struct WideCopy <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
end
SpaceAGORA.bind_spacecraft(::WideCopy, k::Int) = WideCopy(k)
SpaceAGORA.update_control!(::WideCopy, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(::WideCopy, x::StateSample, p, t::Float64) = (Z3, Z3)
struct NarrowCopy <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
end
SpaceAGORA.bind_spacecraft(::NarrowCopy, k::Int) = NarrowCopy(k)
SpaceAGORA.spacecraft_scope(m::NarrowCopy) = m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_control!(::NarrowCopy, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(::NarrowCopy, x::StateSample, p, t::Float64) = (Z3, Z3)

# A control model whose scope is a selection and that has no bind_spacecraft method.
struct SelectsCtl <: SM.AbstractControlEffectorModel end
SpaceAGORA.spacecraft_scope(::SelectsCtl) = SpacecraftSelection(ids=[11])
SpaceAGORA.update_control!(::SelectsCtl, u, p, t::Float64, i::Int) = nothing

# A control model with the scope given, for the ensembles.
struct ScopedCtl <: SM.AbstractControlEffectorModel
    scope::Any
end
SpaceAGORA.spacecraft_scope(m::ScopedCtl) = m.scope
SpaceAGORA.update_control!(::ScopedCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(::ScopedCtl, x::StateSample, p, t::Float64) = (Z3, Z3)

# Two spacecraft for an ensemble, with configuration-level `control` and `declared`
# on each spacecraft.
function ensemble_args(; control=(), declared=(), n=2)
    sats = map(1:n) do k
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(a=PLANET.Rp_e + 700e3, e=0.01, i=45.0, Ω=30.0k, ω=0.0, ν=10.0k)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=10.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=10 + k, control=SM.ControlModel(Tuple(declared), fill(1.0, length(declared))))
    end
    env = SM.EnvironmentModel(planet=PLANET, EI=200.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=5.0, orientation_sim=false, num_steps_to_save=10),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], (SM.InverseSquaredJ2GravityModel(),)),
        control_model=SM.ControlModel(control_effectors=Tuple(control), control_rates=fill(1.0, length(control))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

@testset "Maneuver guidance declared on spacecraft" begin
    @testset "a copy acts on its spacecraft only, and the positional form is kept" begin
        template = raise()
        @test spacecraft_scope(template) === AllSpacecraft()
        bound = SpaceAGORA.bind_spacecraft(template, 2)
        @test bound.spacecraft_index == 2 && spacecraft_scope(bound) == BoundSpacecraft((2,))
        @test bound.command_state !== template.command_state        # its own run state
        @test SpaceAGORA.bind_spacecraft(bound, 2) === bound
        table = campaign([1, 2], [2.0, -1.0])
        copy2 = SpaceAGORA.bind_spacecraft(table, 2)
        @test spacecraft_scope(copy2) == BoundSpacecraft((2,)) && copy2.maneuver_Δv === table.maneuver_Δv
        @test SM.ApoapsisTargetPeriapsisRaiseGuidanceModel(7.0e6, 2.0e5, 0.0, 0.5, Int64[]).spacecraft_index == 0
        @test SM.AerobrakingCampaignPropulsiveManeuverGuidanceModel([1], [1.0], Float64[]).spacecraft_index == 0
    end

    @testset "copies with equal settings give the bits of one model for every spacecraft" begin
        for (scenario, make) in ((:raise, () -> raise()), (:campaign, () -> campaign([1, 2], [2.0, -1.0])))
            @testset "$scenario" begin
                plain = solve(maneuver_args(scenario))
                one = solve(maneuver_args(scenario; shared=(make(),)))
                # Every spacecraft maneuvered: its semimajor axis moved by kilometres.
                @test all(k -> abs(final_orbit(one, k).a - final_orbit(plain, k).a) > 1.0e3, 1:3)
                template = make()
                found = Ref{Any}(nothing)
                logs = warnings(BOUND_ALL) do
                    found[] = solve(maneuver_args(scenario; declared=((template,), (template,), (template,))))
                end
                @test same_bits(one, found[])
                @test isequal(collect(one.u[end]), collect(found[].u[end]))
                @test isempty(logs)                                     # the built-ins' copies are bound
            end
        end
    end

    @testset "copies with their own settings each follow their own" begin
        # Each spacecraft of the mixed run ends where it ends in a run that gives
        # every spacecraft its settings.
        targets = (450e3, 500e3, 550e3)
        mixed = solve(maneuver_args(:raise; declared=map(h -> (raise(h_p=h),), targets)))
        for (k, h) in enumerate(targets)
            alone = solve(maneuver_args(:raise; shared=(raise(h_p=h),)))
            @test isapprox(final_orbit(mixed, k).periapsis, final_orbit(alone, k).periapsis; rtol=1e-6)
        end
        # 50, 100 and 150 km above the 400 km periapsis, to the oblate altitude.
        raised = [final_orbit(mixed, k).periapsis - (PLANET.Rp_e + 400e3) for k in 1:3]
        @test issorted(raised) && all(>(10e3), diff(raised))
        tables = (([1], [2.0]), ([1], [-1.0]), ([1], [1.0]))
        mixed = solve(maneuver_args(:campaign; declared=map(t -> (campaign(t...),), tables)))
        for (k, t) in enumerate(tables)
            alone = solve(maneuver_args(:campaign; shared=(campaign(t...),)))
            @test isapprox(final_orbit(mixed, k).a, final_orbit(alone, k).a; rtol=1e-6)
        end
        plain = solve(maneuver_args(:campaign))
        Δa = [final_orbit(mixed, k).a - final_orbit(plain, k).a for k in 1:3]
        @test Δa[1] > Δa[3] > 1.0e3 && Δa[2] < -1.0e3                     # +2, +1 and -1 m/s at orbit 1
    end

    @testset "bound copies share no run state the run changes" begin
        for model in (raise(), campaign([1], [1.0]))
            fixture = C.conformance_fixture(guidance=(model,))
            check = C._bind_fresh_state_check(model, :guidance_model, C.FIRST_RELEASE;
                configuration=fixture.configuration, index=1)
            @test check.status === :pass
        end
    end

    @testset "a bound copy checkpoints its own entry" begin
        bound = SpaceAGORA.bind_spacecraft(raise(), 2)
        @test L.checkpoint_state(bound) == (entry=nothing,)               # before its first tick
        push!(bound.command_state, 0, 1)
        @test L.checkpoint_state(bound) == (entry=1,)
        fresh = SpaceAGORA.bind_spacecraft(raise(), 2)
        L.restore_checkpoint_state!(fresh, L.checkpoint_state(bound))
        @test fresh.command_state == [0, 1]
        L.restore_checkpoint_state!(fresh, (entry=nothing,))
        @test fresh.command_state == [0, 1]                                # nothing to restore
        shared = raise()
        shared.command_state = [2, 0, 1]
        @test L.checkpoint_state(shared) == (command_state=[2, 0, 1],)     # the form of SpaceAGORA 0.6
    end

    @testset "a resumed run of bound copies is the shared model's resumed run" begin
        # With one template for every spacecraft and with one each, resumed from a
        # checkpoint before the first guidance tick (30 s) and from one after the
        # burns: the bits of the shared model resumed from the same checkpoint, and
        # at the early break those of the run without it.
        checkpointed(dir, interval; resume=false) = (checkpoint_enabled=true, checkpoint_interval_s=interval,
            checkpoint_directory=dir, resume_from_checkpoint=resume)
        function resumed(kw, mission, interval)
            split = mktempdir()
            solve(maneuver_args(; kw..., mission_time=interval, settings=checkpointed(split, interval)))
            return solve(maneuver_args(; kw..., mission_time=mission, settings=checkpointed(split, interval; resume=true)))
        end
        one_template() = (t = raise(); (declared=((t,), (t,), (t,)),))
        own_templates() = (declared=((raise(),), (raise(),), (raise(),)),)
        for (mission, interval) in ((100.0, 10.0), (14000.0, 7000.0))
            shared = resumed((shared=(raise(),),), mission, interval)
            for templates in (one_template, own_templates)
                bound = resumed(templates(), mission, interval)
                @test bound.t[end] == shared.t[end] && isequal(collect(bound.u[end]), collect(shared.u[end]))
            end
            # The uninterrupted run with the same checkpoint interval must match.
            whole = solve(maneuver_args(shared=(raise(),), mission_time=mission,
                settings=checkpointed(mktempdir(), interval)))
            @test isequal(collect(shared.u[end]), collect(whole.u[end]))
        end
        # A campaign run past its first orbit at the 7000 s checkpoint preserves
        # the scheduler's waiting commands and plans and the orbit counter.
        # The #509 checkpoint repair makes the resumed final state exact.
        g = () -> campaign([1, 2], [2.0, -1.0])
        whole = solve(maneuver_args(:campaign; shared=(g(),), settings=checkpointed(mktempdir(), 7000.0)))
        split = mktempdir()
        solve(maneuver_args(:campaign; shared=(g(),), mission_time=7000.0, settings=checkpointed(split, 7000.0)))
        late = solve(maneuver_args(:campaign; shared=(g(),), settings=checkpointed(split, 7000.0; resume=true)))
        @test isequal(collect(late.u[end]), collect(whole.u[end]))
    end

    @testset "the burn scheduler cannot be declared on a spacecraft" begin
        scheduler = SM.ThrusterModels.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0],
            start_burn_time=[-1.0], stop_burn_time=[-1.0], Isp=[300.0])
        @test_throws "one model for every spacecraft with one entry per spacecraft" SpaceAGORA.bind_spacecraft(
            scheduler, 1)
        # The pointer is for a model that acts on every spacecraft, not for a bound one.
        @test_throws "acts on every spacecraft" SpaceAGORA.bind_spacecraft(scheduler, 1)
        # A model whose scope is a selection is told so, without the shared-model pointer.
        message = try
            SpaceAGORA.bind_spacecraft(SelectsCtl(), 1)
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("the spacecraft its spacecraft_scope, a SpacecraftSelection, selects", message)
        @test !occursin("burn scheduler", message) && !occursin("AllSpacecraft", message)
    end

    @testset "a bound copy acts on its spacecraft only, also for a direct caller" begin
        args = maneuver_args(:raise)
        p = SM.ODEParams(n_sats=3, args=args)
        u = SpaceAGORA.SimulationEngine.build_initial_conditions(args)
        bound_raise = SpaceAGORA.bind_spacecraft(raise(), 2)
        for model in (bound_raise, SpaceAGORA.bind_spacecraft(campaign([0], [1.0]), 2))
            SpaceAGORA.update_guidance!(model, u, p, 0.0, 1)
            @test isequal(p.shared_buffers.maneuver_commands[1], SM.PropulsiveManeuverCommand())
        end
        @test isempty(bound_raise.command_state)                          # not grown for spacecraft 1
        # The model for every spacecraft acts on any of them, as before.
        SpaceAGORA.update_guidance!(campaign([0], [1.0]), u, p, 0.0, 1)
        @test p.shared_buffers.maneuver_commands[1].valid
    end

    @testset "a bound copy that acts on every spacecraft warns once per type" begin
        wide = warnings(BOUND_ALL) do
            run_simulation(ensemble_args(declared=(WideCopy(0),)))
        end
        @test length(wide) == 1 &&
            occursin("WideCopy is declared on the spacecraft at run position 1 (id 11)", string(only(wide).message))
        @test isempty(warnings(() -> run_simulation(ensemble_args(declared=(NarrowCopy(0),))), BOUND_ALL))
    end

    @testset "an ensemble refuses configuration-level GNC bound to run positions" begin
        message = "each ensemble member runs its spacecraft alone"
        for scope in (BoundSpacecraft((1,)), CoupledSpacecraft(1, (2,)), SpacecraftSelection(ids=[12]),
                SpacecraftSelection(positions=[1]))
            @test_throws message run_constellation_ensemble(ensemble_args(control=(ScopedCtl(scope),));
                allow_gnc_effectors=true)
        end
        @test_throws "control model 1 (ScopedCtl)" run_constellation_ensemble(
            ensemble_args(control=(ScopedCtl(BoundSpacecraft((1,))),)); allow_gnc_effectors=true)
        # A model for every spacecraft, and models declared on the spacecraft, run.
        all_members = run_constellation_ensemble(ensemble_args(control=(ScopedCtl(AllSpacecraft()),));
            allow_gnc_effectors=true)
        @test length(all_members.successful) == 2
        declared = run_constellation_ensemble(ensemble_args(declared=(NarrowCopy(0),)))
        @test length(declared.successful) == 2
        # An ensemble of one keeps a scope that names position 1 alone, as a run would.
        for scope in (BoundSpacecraft((1,)), SpacecraftSelection(ids=[11]))
            @test length(run_constellation_ensemble(ensemble_args(control=(ScopedCtl(scope),), n=1);
                allow_gnc_effectors=true).successful) == 1
        end
        @test_throws message run_constellation_ensemble(ensemble_args(control=(ScopedCtl(BoundSpacecraft((2,))),), n=1);
            allow_gnc_effectors=true)
        # The scope as written, without its type parameters.
        @test_throws "declares spacecraft_scope BoundSpacecraft((1,))," run_constellation_ensemble(
            ensemble_args(control=(ScopedCtl(BoundSpacecraft((1,))),)); allow_gnc_effectors=true)
    end
end

end
