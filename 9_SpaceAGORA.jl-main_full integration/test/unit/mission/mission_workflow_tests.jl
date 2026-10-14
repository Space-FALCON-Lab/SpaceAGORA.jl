module MissionWorkflowTests

using Test
using SpaceAGORA
using OrdinaryDiffEq: DiscreteCallback
using DataFrames
using JSON
const SM = SpaceAGORA.SimulationModel
const MW = SpaceAGORA.MissionWorkflow

function configuration(; duration_s=30.0)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
        prop_mass=0.0, id=1)
    return make_example_config(; planet, spacecraft, mission_time=Float64(duration_s),
        initial_time=InitialTime(year=2024, month=1, day=1),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, verbose=false, results=false,
        solver_config=SolverConfig(solver_mode=:tsit5, parallel=false))
end

function with_settings(config; kwargs...)
    settings = config.simulation_settings
    fields = fieldnames(typeof(settings))
    values = NamedTuple{fields}(map(n -> getfield(settings, n), fields))
    return SM.SimConfig._with_configuration(config;
        simulation_settings=SimulationSettings(; merge(values, (; kwargs...))...))
end

function assert_same_table(a, b)
    @test names(a) == names(b)
    @test size(a) == size(b)
    for name in names(a)
        @test isequal(a[!, name], b[!, name])
    end
end

@testset "Mission definitions are bounded and descriptive" begin
    @test_throws ArgumentError MissionPhase(:bad; window_s=(-1.0, 1.0))
    @test_throws ArgumentError MissionPhase(:bad; window_s=(1.0, 1.0))
    @test_throws ArgumentError MissionPhase(:bad; window_s=(2.0, 1.0))
    @test_throws ArgumentError MissionPhase(:bad; window_s=(0.0, Inf))
    @test_throws ArgumentError MissionPhase(:bad; window_s=(NaN, 1.0))
    @test_throws ArgumentError MissionPhase(:bad; window_s=(0.0,))
    @test_throws ArgumentError MissionPhase(" "; window_s=(0.0, 1.0))
    phase = MissionPhase(:coast; window_s=(0, 30))
    @test phase.window_s === (0.0, 30.0)
    @test_throws ArgumentError MissionPlan(" ", configuration)
    @test_throws ArgumentError MissionPlan("bad", configuration; phases=(nothing,))
    @test_throws ArgumentError MissionPlan("bad", configuration; events=(nothing,))
    @test_throws ArgumentError MissionPlan("bad", configuration; phases=(phase, phase))
    @test_throws ArgumentError MissionPlan("bad", configuration; phases=(phase, MissionPhase(:late; window_s=(20, 40))))
    @test_throws ArgumentError MissionPlan("bad", configuration; phases=(MissionPhase(:late; window_s=(40, 50)), phase))
    event = MissionEvent(:sample, _ -> DiscreteCallback((u,t,i) -> false, _ -> nothing))
    @test_throws ArgumentError MissionPlan("bad", configuration; events=(event, event))
    @test_throws ArgumentError MissionEvent(" ", identity)
    plan = MissionPlan("coast", configuration; phases=[phase], events=[event])
    @test plan.phases === (phase,)
    @test plan.events === (event,)
end

@testset "Mission validation creates no outputs" begin
    mktempdir() do root
        destination = joinpath(root, "not-created")
        plan = MissionPlan("coast", configuration)
        @test_throws ArgumentError run_mission(plan; output_dir=destination, save_scene=false)
        @test !ispath(destination)
        @test_throws ArgumentError run_mission(MissionPlan("wrong", () -> nothing); output_dir=destination)
        @test !ispath(destination)
        @test_throws ErrorException run_mission(MissionPlan("throws", () -> error("builder failed")); output_dir=destination)
        @test !ispath(destination)
        too_long = MissionPlan("long", configuration; phases=(MissionPhase(:beyond; window_s=(0, 31)),))
        @test_throws ArgumentError run_mission(too_long; output_dir=destination)
        @test !ispath(destination)
        for setting in (:checkpoint_enabled, :resume_from_checkpoint)
            checkpointed = MissionPlan("checkpoint", () -> with_settings(configuration(); NamedTuple{(setting,)}((true,))...))
            @test_throws ArgumentError run_mission(checkpointed; output_dir=destination)
            @test !ispath(destination)
        end
        withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
            @test_throws ArgumentError run_mission(plan; output_dir=destination)
            @test !ispath(destination)
        end
        @test_throws ArgumentError run_mission(plan; output_dir=" ")
        occupied = joinpath(root, "occupied")
        mkdir(occupied)
        sentinel = joinpath(occupied, "keep.txt")
        write(sentinel, "keep exactly")
        @test_throws ArgumentError run_mission(plan; output_dir=occupied)
        @test read(sentinel, String) == "keep exactly"
        @test readdir(occupied) == ["keep.txt"]
        @test_throws ArgumentError run_mission(plan; output_dir=sentinel)
        alias = joinpath(root, "alias")
        symlink(occupied, alias)
        @test_throws ArgumentError run_mission(plan; output_dir=alias)
        for (name, factory) in ((:bad_value, _ -> 42), (:factory_failure, _ -> error("factory failed")),
                                 (:wrong_arity, () -> nothing))
            bad = MissionPlan("bad factory", configuration; events=(MissionEvent(name, factory),))
            exception = name == :factory_failure ? ErrorException : ArgumentError
            @test_throws exception run_mission(bad; output_dir=destination, save_scene=false, save_payload=false)
            @test !ispath(destination)
        end
    end
end

@testset "Mission output choices preserve all scientific configuration fields" begin
    built = Ref{Any}()
    plan = MissionPlan("inspection", () -> (built[] = configuration()))
    prepared = MW._configuration(plan, (;), "/unused/mission-output", true)
    original = built[]
    for field in fieldnames(typeof(original))
        field == :simulation_settings && continue
        @test getfield(prepared, field) === getfield(original, field)
    end
    for field in fieldnames(typeof(original.simulation_settings))
        field in (:results, :results_directory, :save_visualization_scene) && continue
        @test getfield(prepared.simulation_settings, field) === getfield(original.simulation_settings, field)
    end
    @test prepared.simulation_settings.results
    @test prepared.simulation_settings.save_visualization_scene
    @test prepared.simulation_settings.results_directory == "/unused/mission-output"
    @test !original.simulation_settings.results
end

@testset "One solve: direct and mission numerical parity" begin
    mktempdir() do root
        direct_config = with_results_directory(configuration(), joinpath(root, "direct"); visualization=false)
        direct = run_simulation(direct_config; return_results=true, return_solution=true, visualization=false)
        builder_calls = Ref(0)
        builder = (; duration_s=30.0) -> begin
            builder_calls[] += 1
            configuration(; duration_s)
        end
        plan = MissionPlan("coast", builder;
            phases=(MissionPhase(:early; window_s=(0, 13)), MissionPhase(:late; window_s=(13, 30))))
        mission = run_mission(plan; output_dir=joinpath(root, "mission"), parameters=(duration_s=30.0,),
            save_scene=false, save_payload=false, return_solution=true)
        @test builder_calls[] == 1
        @test mission isa MissionRun
        @test mission.name == "coast"
        @test mission.payload_path === nothing
        @test mission.results.solution.t == direct.solution.t # descriptive windows insert no stops
        @test isequal(mission.results.solution.u, direct.solution.u)
        assert_same_table(mission.results.table, direct.table)
        @test mission.phases === plan.phases
        @test isempty(mission.event_names)
        @test !isempty(mission.files)
        @test all(isfile, mission.files)
        @test all(path -> dirname(path) == mission.directory, mission.files)
        @test mission.files == mission.results.files
        @test_throws ArgumentError run_mission(plan; output_dir=mission.directory, save_scene=false, save_payload=false)
    end
end

@testset "Event factories see isolated configuration and own fresh callback state" begin
    mktempdir() do root
        built = Any[]
        captured = Any[]
        counters = Base.RefValue{Int}[]
        fired = Tuple{Int, Float64, Bool}[]
        builder = () -> begin
            config = configuration()
            push!(built, config)
            config
        end
        factory = config -> begin
            push!(captured, config)
            count = Ref(0)
            push!(counters, count)
            call = length(counters)
            condition(u, t, integrator) = t >= 5.0 && count[] == 0
            affect!(integrator) = begin
                count[] += 1
                push!(fired, (call, Float64(integrator.t), config === integrator.p.args))
                nothing
            end
            DiscreteCallback(condition, affect!; save_positions=(false, false))
        end
        plan = MissionPlan("event", builder; events=(MissionEvent(:observe, factory),))
        a = run_mission(plan; output_dir=joinpath(root, "one"), save_scene=false, save_payload=false)
        b = run_mission(plan; output_dir=joinpath(root, "two"), save_scene=false, save_payload=false)
        @test length(built) == length(captured) == length(counters) == 2
        @test built[1] !== built[2]
        @test captured[1] === a.results.configuration
        @test captured[2] === b.results.configuration
        @test captured[1].dynamics_model.spacecraft !== built[1].dynamics_model.spacecraft
        @test captured[2].dynamics_model.spacecraft !== built[2].dynamics_model.spacecraft
        @test counters[1] !== counters[2]
        @test getindex.(counters) == [1, 1]
        @test length(fired) == 2
        @test all(last, fired)
        @test fired[1][2] == fired[2][2]
        @test a.event_names == (:observe,)
        assert_same_table(a.results.table, b.results.table)
    end
end

@testset "Existing callbacks precede factory callbacks" begin
    mktempdir() do root
        order = Symbol[]
        one_shot = name -> begin
            seen = Ref(false)
            DiscreteCallback((u,t,i) -> t >= 5.0 && !seen[], i -> begin
                seen[] = true
                push!(order, name)
            end; save_positions=(false, false))
        end
        config = with_results_directory(configuration(), joinpath(root, "ordered"); visualization=false)
        run_simulation(config; extra_callbacks=(one_shot(:existing),),
            extra_callback_factories=(_ -> one_shot(:factory),), visualization=false)
        @test order == [:existing, :factory]
    end
end

@testset "Saved mission payload can be moved and viewed without rerunning" begin
    mktempdir() do root
        calls = Ref(0)
        plan = MissionPlan("portable", () -> begin
            calls[] += 1
            configuration()
        end)
        run = run_mission(plan; output_dir=joinpath(root, "original"))
        @test calls[] == 1
        @test run.payload_path !== nothing
        @test isfile(run.payload_path)
        @test run.payload_path in run.files
        @test isfile(run.prefix * "_scene.json")
        @test isfile(run.prefix * ".feather")
        payload = read(run.payload_path)
        @test JSON.parse(String(copy(payload))) isa AbstractDict
        moved = joinpath(root, "moved")
        cp(run.directory, moved)
        rm(run.directory; recursive=true)
        relocated_payload = joinpath(moved, basename(run.payload_path))
        @test read(relocated_payload) == payload
        # Rebuild from the relocated saved scene + Arrow bundle without a live configuration.
        exported = SM.SceneVisualization.export_viewer_payload(joinpath(moved, basename(run.prefix));
            out=joinpath(moved, "second_viewer.json"), textures=false)
        @test isfile(exported)
        @test JSON.parsefile(exported) == JSON.parsefile(relocated_payload)
        @test calls[] == 1
        # Renderer availability is checked before reading the old handle path.
        # In a full integration process the renderer is already loaded; a
        # core-only process instead reports the actionable optional-package error.
        expected_error = visualizer_available() ? ArgumentError : VisualizerUnavailableError
        @test_throws expected_error export_visualization(run; textures=false)
    end
end

end # module MissionWorkflowTests
