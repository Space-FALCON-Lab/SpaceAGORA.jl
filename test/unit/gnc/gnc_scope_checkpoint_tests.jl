# The GNC models' spacecraft scopes in checkpoint files (`gnc_scopes`): a resumed
# run compares the scopes its guidance, navigation and control models have after
# the run-start initializers with the checkpoint's, before any model's run state is
# restored. A model bound to the same spacecraft resumes as an uninterrupted run, bit
# for bit; one bound to another spacecraft (an index-bound model, a bound copy, a
# selection) is refused, with or without run state of its own. A run with
# guidance, navigation or control models refuses a checkpoint written without the
# entry, by SpaceAGORA 0.5.1 or earlier, before anything is restored; a run
# without them resumes from it. A damaged entry is refused as it is read, and a
# bound scope that names positions outside the run is recorded as its ticks
# visit it, so checkpointing such a run changes nothing.
module GNCScopeCheckpointTests

using Test, StaticArrays, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const C = SpaceAGORA.Conformance
const ZERO3 = SVector(0.0, 0.0, 0.0)

# A momentum manager whose command turns its spacecraft, as in the built-in models'
# checkpoint tests: it keeps run state (the wheel momentum and the held torque).
manager(index=1) = SM.ControlHooks.MagneticMomentumManagerModel(spacecraft_index=index,
    h_wheels_0=SVector(0.05, -0.02, 0.01), commanded_torque=(t, r, v, q, w) -> SVector(1.0e-5, -2.0e-5, 5.0e-6),
    b_field_ii=(t, r) -> SVector(2.0e-5, 1.0e-5, -3.0e-5))

# A model bound to one spacecraft, pushing it, that keeps no run state
# (checkpoint_state is nothing).
struct BoundPush <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
end
SpaceAGORA.update_control!(::BoundPush, u, p, t, i) = nothing
SpaceAGORA.control_wrench(m::BoundPush, x::StateSample, p, t) =
    x.spacecraft_index == m.spacecraft_index ? (SVector(0.0, 1.0e-3, 0.0), ZERO3) : (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::BoundPush) = BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.bind_spacecraft(::BoundPush, i::Int) = BoundPush(i)

# A model bound to spacecraft 2 whose initializer widens its scope past the run,
# after the pre-run check: its ticks, and its push, reach spacecraft 2 only.
mutable struct Widened <: SM.AbstractControlEffectorModel
    indices::Tuple{Vararg{Int}}
end
SpaceAGORA.update_control!(::Widened, u, p, t, i) = nothing
SpaceAGORA.control_wrench(m::Widened, x::StateSample, p, t) =
    x.spacecraft_index in m.indices ? (SVector(0.0, 1.0e-3, 0.0), ZERO3) : (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::Widened) = BoundSpacecraft(m.indices)
SpaceAGORA.initialize_control!(m::Widened, u, p, t) = (m.indices = (2, 99); nothing)

# Models acting on every spacecraft that keep no run state.
struct Idle <: SM.AbstractControlEffectorModel end
SpaceAGORA.update_control!(::Idle, u, p, t, i) = nothing
SpaceAGORA.control_wrench(::Idle, x::StateSample, p, t) = (ZERO3, ZERO3)
struct IdleGuidance <: AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::IdleGuidance, u, p, t::Float64, i::Int) = nothing
struct IdleNavigation <: AbstractNavigationModel end
SpaceAGORA.update_navigation!(::IdleNavigation, u, p, t::Float64, i::Int) = nothing

# A model with whatever scope it is given (guidance, navigation and control models
# do not take a SpacecraftSelection in a run yet; the scopes below are compared
# without one).
struct Selected <: SM.AbstractControlEffectorModel
    scope::Any
end
SpaceAGORA.spacecraft_scope(m::Selected) = m.scope

# Two spacecraft on the conformance fixture's orbit, with attitude, each model
# ticking every second: `control`, `guidance` and `navigation` at configuration
# level, and the control models `declared[i]` on spacecraft `i` (bound to it when
# the run starts).
function scope_args(; control=(), guidance=(), navigation=(), declared=((), ()), mission_time=6.0, settings=(;))
    planet = SM.make_no_gram_planet(:earth)
    spacecraft = map(1:2) do id
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        orbit = SM.InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 650e3, i=45.0, ω=0.0, Ω=0.0,
            ν=10.0 + 5.0 * (id - 1))
        mine = Tuple(declared[id])
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
            initial_condition=orbit, id=id, control=SM.ControlModel(mine, fill(1.0, length(mine))))
    end
    environment = SM.EnvironmentModel(planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(; results=false, verbose=false, generate_plots=false, normalize=false,
            settings...),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=true, num_steps_to_save=100),
        environment_model=environment,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[spacecraft...], (SM.InverseSquaredJ2GravityModel(),)),
        guidance_model=SM.GuidanceModel(Tuple(guidance), fill(1.0, length(guidance))),
        navigation_model=SM.NavigationModel(Tuple(navigation), fill(1.0, length(navigation))),
        control_model=SM.ControlModel(Tuple(control), fill(1.0, length(control))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

checkpointed(dir; resume=false) =
    (checkpoint_enabled=true, checkpoint_interval_s=1.5, checkpoint_directory=dir, resume_from_checkpoint=resume)
# A run of 3 s that writes its checkpoints into `dir`, and one of 6 s resumed from them.
first_half(dir; kw...) = scope_args(; mission_time=3.0, settings=checkpointed(dir), kw...)
second_half(dir; kw...) = scope_args(; mission_time=6.0, settings=checkpointed(dir; resume=true), kw...)

# The message of the `ArgumentError` that `f()` refuses with, or `nothing` when it returns.
function refusal(f)
    try
        f()
        return nothing
    catch e
        e isa ArgumentError || rethrow()
        return sprint(showerror, e)
    end
end

saved_runtime(args) = SM.IOSerialization._load_checkpoint(args).runtime_state

# Rewrites the checkpoint in `args`' directory with `change(runtime_state)` as its
# run state. Dropping `gnc_scopes` gives the file SpaceAGORA 0.5.1 wrote for the run.
function rewrite_checkpoint!(change, args)
    ck = SM.IOSerialization._load_checkpoint(args)
    SM.IOSerialization._write_checkpoint!(args, ck.t, ck.u, SE.CHECKPOINT_SCHEMA_VERSION;
        solver_mode=ck.solver_mode, runtime_state=change(ck.runtime_state))
    return nothing
end
without_scopes(runtime) = Base.structdiff(runtime, NamedTuple{(:gnc_scopes,)})
# A copy of the checkpoints in `dir` as SpaceAGORA 0.5.1 wrote them. A run resumed
# from a directory writes its own checkpoints there, so each resume takes a copy.
function legacy_copy(dir)
    legacy = mktempdir()
    cp(dir, legacy; force=true)
    rewrite_checkpoint!(without_scopes, first_half(legacy))
    return legacy
end

# The kit's checkpoint round trip, which compares a run resumed from its checkpoint
# halfway with an uninterrupted run with the same checkpoints, bit for bit: the
# saved states and the control model at `index` of the run's flattened models.
round_trip(configuration, subject; index=1) = C._checkpoint_round_trip_check(
    (cfg, factories) -> C._model_solve(cfg, :control_model, index; factories), C._model_run_difference, subject,
    configuration)

@testset "GNC spacecraft scopes in checkpoint files" begin
    @testset "a scope is plain data taken afresh, the same for the same spacecraft" begin
        args = scope_args()
        with_control(models...) = SM.SimConfig._with_configuration(args;
            control_model=SM.ControlModel(Tuple(models), fill(1.0, length(models))))
        scopes(models...) = SE._gnc_scope_baseline(with_control(models...)).control
        @test SE._gnc_scope_baseline(args) == (version=1, navigation=Any[], guidance=Any[], control=Any[])
        @test scopes(Selected(AllSpacecraft())) == [(scope=:all,)]
        @test scopes(Selected(CoupledSpacecraft(2, (1,)))) == [(scope=:coupled, primary=2, partners=[1])]
        # A bound scope is the run positions it names, sorted and without repeats,
        # however it names them; naming every spacecraft is still bound.
        bound = [(scope=:bound, positions=[1])]
        for same in (BoundSpacecraft((1,)), BoundSpacecraft((1, 1)), SpacecraftSelection(positions=1),
                     SpacecraftSelection(ids=1), SpacecraftSelection(where=(i, sc) -> sc.id == 1))
            @test scopes(Selected(same)) == bound
        end
        @test scopes(Selected(BoundSpacecraft((2, 1)))) == [(scope=:bound, positions=[1, 2])]
        # Positions outside the run's two are never ticked, and not recorded.
        @test scopes(Selected(BoundSpacecraft((3, 2, 0, 2)))) == [(scope=:bound, positions=[2])]
        @test scopes(Selected(BoundSpacecraft((3,)))) == [(scope=:bound, positions=Int[])]
        @test !isequal(only(scopes(Selected(BoundSpacecraft((1, 2))))), only(scopes(Selected(AllSpacecraft()))))
        # A selection changed in place still equals itself, as the selection a model
        # holds would; the baseline taken before keeps the spacecraft it named, so
        # the two compare unequal and the resume is refused.
        selection = SpacecraftSelection(positions=[1])
        model = Selected(selection)
        taken = SE._gnc_scope_baseline(with_control(model))
        held = model.scope
        selection.positions[1] = 2
        @test held == selection
        @test taken.control == bound
        now = SE._gnc_scope_baseline(with_control(model))
        @test now.control == [(scope=:bound, positions=[2])]
        msg = refusal(() -> SE._check_gnc_scope_baseline(with_control(model), taken, now))
        @test msg !== nothing && occursin("control model 1 (Selected) acts on the spacecraft at run position 2 here " *
            "and on the spacecraft at run position 1 in the checkpoint", msg)
        @test SE._check_gnc_scope_baseline(with_control(model), now, now) === nothing
        # The other forms name their spacecraft in the refusal.
        coupled = SE._gnc_scope_baseline(with_control(Selected(CoupledSpacecraft(2, (1,)))))
        msg = refusal(() -> SE._check_gnc_scope_baseline(with_control(Selected(AllSpacecraft())), coupled,
            SE._gnc_scope_baseline(with_control(Selected(AllSpacecraft())))))
        @test msg !== nothing && occursin("acts on every spacecraft here and on the spacecraft at run position 2, " *
            "coupled with run position 1 in the checkpoint", msg)
    end

    @testset "the same scopes resume as an uninterrupted run, bit for bit" begin
        # The manager of spacecraft 2 at configuration level, the same manager
        # declared on spacecraft 2 (a bound copy), and a bound model without run
        # state of its own.
        configured = scope_args(control=(manager(2),))
        declared = scope_args(declared=((), (manager(1),)))
        stateless = scope_args(control=(BoundPush(2),))
        for (name, configuration, subject) in (("configured", configured, manager(2)),
                                              ("declared", declared, manager(1)), ("stateless", stateless, BoundPush(2)))
            check = round_trip(configuration, subject)
            @testset let case = name, status = check.status, detail = check.detail
                @test status === :pass
            end
        end
        # The checkpoint records the scope the run used, the same for both declarations.
        for (configuration, kept) in ((configured, true), (declared, true), (stateless, false))
            args = SM.SimConfig._with_configuration(configuration; simulation_settings=SM.SimulationSettings(
                results=false, verbose=false, generate_plots=false, normalize=false, checkpoint_enabled=true,
                checkpoint_interval_s=3.0, checkpoint_directory=mktempdir()))
            run_simulation(args)
            runtime = saved_runtime(args)
            @test runtime.gnc_scopes ==
                (version=1, navigation=Any[], guidance=Any[], control=Any[(scope=:bound, positions=[2])])
            @test hasproperty(runtime, :models) == kept             # nothing kept: the scopes alone
        end
    end

    @testset "a model acting on other spacecraft is refused before any run state is restored" begin
        # The checkpoint holds the run state of the manager of spacecraft 1.
        dir = mktempdir()
        written = first_half(dir; control=(manager(1),))
        run_simulation(written)
        saved = only(saved_runtime(written).models.control)
        @test saved.ticks > 0 && saved.held_torque_body != ZERO3
        # Resumed with the manager of spacecraft 2: refused after its initializer,
        # which reset it, and before its run state came back.
        moved = manager(2)
        msg = refusal(() -> run_simulation(second_half(dir; control=(moved,)); isolate_state=false))
        @test msg !== nothing && occursin("control model 1 (MagneticMomentumManagerModel) acts on the spacecraft at " *
            "run position 2 here and on the spacecraft at run position 1 in the checkpoint", msg)
        @test moved.ticks == 0 && moved.held_torque_body == ZERO3 && moved.run_written === nothing
        # A bound copy: the manager declared on spacecraft 1, then on spacecraft 2.
        dir = mktempdir()
        run_simulation(first_half(dir; declared=((manager(),), ())))
        msg = refusal(() -> run_simulation(second_half(dir; declared=((), (manager(),)))))
        @test msg !== nothing && occursin("control model 1 (MagneticMomentumManagerModel) acts on the spacecraft at " *
            "run position 2 here and on the spacecraft at run position 1 in the checkpoint", msg)
        @test refusal(() -> run_simulation(second_half(dir; declared=((manager(),), ())))) === nothing
        # A model without run state, so with no `models` entry to compare.
        dir = mktempdir()
        written = first_half(dir; control=(BoundPush(1),))
        run_simulation(written)
        @test !hasproperty(saved_runtime(written), :models)
        msg = refusal(() -> run_simulation(second_half(dir; control=(BoundPush(2),))))
        @test msg !== nothing && occursin("control model 1 (BoundPush) acts on the spacecraft at run position 2 here " *
            "and on the spacecraft at run position 1 in the checkpoint", msg)
        # Another number of models, whose scopes cannot be paired.
        msg = refusal(() -> run_simulation(second_half(dir; control=(BoundPush(1), BoundPush(2)))))
        @test msg !== nothing && occursin("it has 1 control model and this run 2", msg)
        @test refusal(() -> run_simulation(second_half(dir; control=(BoundPush(1),)))) === nothing
        # A run without GNC models records that, so its checkpoint is not mistaken
        # for one written before the entry.
        dir = mktempdir()
        written = first_half(dir)
        run_simulation(written)
        @test saved_runtime(written).gnc_scopes == (version=1, navigation=Any[], guidance=Any[], control=Any[])
        msg = refusal(() -> run_simulation(second_half(dir; control=(BoundPush(1),))))
        @test msg !== nothing && occursin("it has 0 control models and this run 1", msg)
        # An entry of another version.
        rewrite_checkpoint!(runtime -> merge(runtime, (gnc_scopes=merge(runtime.gnc_scopes, (version=2,)),)), written)
        msg = refusal(() -> run_simulation(second_half(dir)))
        @test msg !== nothing && occursin("GNC scope checkpoint state version 1 is required", msg)
        # A damaged scope in an entry of the right version is refused when the
        # checkpoint is read, as an ArgumentError naming it.
        dir = mktempdir()
        run_simulation(first_half(dir; control=(BoundPush(1),)))
        for damaged in ((kind=:bound,), (scope=:bound,), (scope=:bound, positions=[1.0]), (scope=:coupled, primary=1),
                        (scope=:coupled, primary=1, partners=nothing), (scope=:other,), :all)
            damaged_dir = mktempdir()
            cp(dir, damaged_dir; force=true)
            rewrite_checkpoint!(runtime -> merge(runtime, (gnc_scopes=merge(runtime.gnc_scopes,
                (control=Any[damaged],)),)), first_half(damaged_dir))
            msg = refusal(() -> run_simulation(second_half(damaged_dir; control=(BoundPush(1),))))
            @testset let damaged = damaged
                @test msg !== nothing && occursin("The checkpoint's GNC scope entry is damaged: it gives control " *
                    "model 1 the scope $(repr(damaged))", msg)
            end
        end
        @test refusal(() -> run_simulation(second_half(dir; control=(BoundPush(1),)))) === nothing
    end

    @testset "a bound scope an initializer widens past the run is recorded as its ticks see it" begin
        # The pre-run check sees spacecraft 2 only; the initializer then names 99 too.
        # The checkpointed run starts, records spacecraft 2, and resumes as an
        # uninterrupted run, bit for bit, as a run without checkpoints starts.
        @test refusal(() -> run_simulation(scope_args(control=(Widened((2,)),)))) === nothing
        check = round_trip(scope_args(control=(Widened((2,)),)), Widened((2,)))
        @testset let status = check.status, detail = check.detail
            @test status === :pass
        end
        dir = mktempdir()
        written = first_half(dir; control=(Widened((2,)),))
        run_simulation(written)
        @test saved_runtime(written).gnc_scopes.control == Any[(scope=:bound, positions=[2])]
        @test refusal(() -> run_simulation(second_half(dir; control=(Widened((2,)),)))) === nothing
    end

    @testset "the scopes are compared as the run-start initializers leave them" begin
        # The RPO planner and controller set their run positions from the spacecraft
        # ids at the run's start. Given the other positions before the run, they are
        # recorded, and compared on a resume, as resolved.
        function rpo(mission_time_s, dir; resume=false)
            args = make_rpo_configuration(planner=DirectRPOPlanner(
                    segment_clearance_at=SpaceAGORA.RPOPlannerLifecycle.station_segment_clearance),
                mission_time_s=mission_time_s, simulation_settings=SM.SimulationSettings(results=false, verbose=false,
                    generate_plots=false, checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir,
                    resume_from_checkpoint=resume))
            for model in (only(args.guidance_model.guidance_effectors), only(args.control_model.control_effectors))
                model.chaser_index, model.target_index = 2, 1
            end
            return args
        end
        dir = mktempdir()
        written = rpo(2.0, dir)
        run_simulation(written)
        resolved = Any[(scope=:coupled, primary=1, partners=[2])]
        @test saved_runtime(written).gnc_scopes == (version=1, navigation=Any[], guidance=resolved, control=resolved)
        @test refusal(() -> run_simulation(rpo(4.0, dir; resume=true))) === nothing
    end

    @testset "a run with GNC models refuses a checkpoint written without the entry" begin
        # The checkpoint of the manager of spacecraft 2, as written now and as
        # SpaceAGORA 0.5.1 wrote it.
        dir = mktempdir()
        run_simulation(first_half(dir; control=(manager(2),)))
        legacy = legacy_copy(dir)
        @test !hasproperty(saved_runtime(first_half(legacy)), :gnc_scopes)
        @test hasproperty(saved_runtime(first_half(legacy)), :models)
        @test refusal(() -> run_simulation(second_half(dir; control=(manager(2),)))) === nothing
        # Refused when the checkpoint is read, before anything is restored, whatever
        # the models: the same manager, another spacecraft's, a stateless bound model,
        # one acting on every spacecraft, guidance or navigation alone. The run's
        # models are left as they were, not even reset by their initializers.
        for (name, models) in (("same", (control=(manager(2),),)), ("moved", (control=(manager(1),),)),
                               ("stateless", (control=(BoundPush(2),),)), ("every spacecraft", (control=(Idle(),),)),
                               ("guidance", (guidance=(IdleGuidance(),),)),
                               ("navigation", (navigation=(IdleNavigation(),),)))
            args = second_half(legacy; models...)
            before = deepcopy(args)
            msg = refusal(() -> run_simulation(args; isolate_state=false))
            @testset let case = name
                @test msg !== nothing && occursin("The checkpoint records no spacecraft scopes for guidance, " *
                    "navigation and control models: it was written by SpaceAGORA 0.5.1 or earlier", msg) &&
                    occursin("Resume it with SpaceAGORA 0.5.1, or run from the start.", msg)
                @test C._same_state(args.navigation_model, before.navigation_model) &&
                    C._same_state(args.guidance_model, before.guidance_model) &&
                    C._same_state(args.control_model, before.control_model)
            end
        end
        moved = manager(2)
        msg = refusal(() -> run_simulation(second_half(legacy; control=(moved,)); isolate_state=false))
        @test occursin("so this run's control model 1 (MagneticMomentumManagerModel) cannot be checked", msg)
        @test moved.ticks == 0 && moved.run_start === nothing                  # its initializer did not run
        # More than three models are named three at a time.
        msg = refusal(() -> run_simulation(second_half(legacy; navigation=(IdleNavigation(),),
            guidance=(IdleGuidance(),), control=(manager(2), BoundPush(1)))))
        @test occursin("navigation model 1 (IdleNavigation), guidance model 1 (IdleGuidance), control model 1 " *
            "(MagneticMomentumManagerModel) and 1 more cannot be checked", msg)
        # A run without guidance, navigation or control models has nothing to
        # compare: it resumes from such a checkpoint as from one written now.
        dir = mktempdir()
        run_simulation(first_half(dir))
        final(args) = run_simulation(args; return_solution=true).u[end]
        from_legacy = final(second_half(legacy_copy(dir)))
        @test isequal(from_legacy, final(second_half(dir)))
    end
end

end
