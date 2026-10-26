module OutputPolicyTests

using Test
using SpaceAGORA
using DataFrames
using Arrow
using CSV
using TOML
using SHA

const SM = SpaceAGORA.SimulationModel

function configuration(directory; results=true, save_csv=true, save_scene=false)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
        prop_mass=0.0, id=1)
    config = make_example_config(; planet, spacecraft, mission_time=2.0,
        initial_time=InitialTime(year=2024, month=1, day=1),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, verbose=false, results,
        results_directory=directory, solver_config=SolverConfig(solver_mode=:tsit5, parallel=false))
    settings = config.simulation_settings
    fields = fieldnames(typeof(settings))
    values = NamedTuple{fields}(map(n -> getfield(settings, n), fields))
    return SM.SimConfig._with_configuration(config;
        simulation_settings=SimulationSettings(; merge(values, (; save_csv, save_visualization_scene=save_scene))...))
end

@testset "Output switches govern files without changing the returned table" begin
    mktempdir() do root
        reference = nothing
        for results in (true, false), save_csv in (true, false), bundle in (true, false)
            directory = joinpath(root, "$results-$save_csv-$bundle")
            mkpath(directory)
            result = withenv("SPACEAGORA_SAVE_BUNDLE" => bundle ? "1" : "0") do
                run_simulation(configuration(directory; results, save_csv); return_results=true)
            end
            @test nrow(result.table) > 1
            if reference === nothing
                reference = result.table
            else
                @test isequal(result.table, reference)
            end
            prefix = joinpath(directory, "simulation_results")
            expected = String[]
            if results && save_csv
                push!(expected, prefix * ".csv")
                @test isequal(CSV.read(prefix * ".csv", DataFrame), result.table)
            end
            if results && bundle
                append!(expected, [prefix * ".feather", prefix * ".manifest.toml"])
                @test isequal(DataFrame(Arrow.Table(prefix * ".feather")), result.table)
                manifest = TOML.parsefile(prefix * ".manifest.toml")
                @test Set(keys(manifest["files"])) == (save_csv ? Set(["csv", "feather"]) : Set(["feather"]))
                for file in values(manifest["files"])
                    @test file["sha256"] == bytes2hex(sha256(read(file["path"])))
                    @test file["size_bytes"] == filesize(file["path"])
                end
            end
            @test Set(result.files) == Set(expected)
            @test Set(readdir(directory; join=true)) == Set(expected)
        end
    end
end

@testset "Disabling CSV preserves old files without claiming them as new output" begin
    mktempdir() do directory
        csv = joinpath(directory, "simulation_results.csv")
        write(csv, "previous run: keep exactly\n")
        before = read(csv)
        for bundle in (true, false)
            result = withenv("SPACEAGORA_SAVE_BUNDLE" => bundle ? "1" : "0") do
                run_simulation(configuration(directory; save_csv=false); return_results=true)
            end
            @test read(csv) == before
            @test csv ∉ result.files
            @test all(!endswith(path, ".csv") for path in result.files)
            expected = bundle ? [joinpath(directory, "simulation_results.feather"),
                joinpath(directory, "simulation_results.manifest.toml")] : String[]
            @test Set(result.files) == Set(expected)
            @test !haskey(TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))["files"], "csv")
        end
        @test withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
            run_simulation(configuration(directory; save_csv=false)) === nothing
        end
    end
end

@testset "Feather-only replacement invalidates previous mission ownership" begin
    mktempdir() do directory
        config = configuration(directory; save_csv=false)
        prefix = SM.IOConfig._results_bundle_prefix(config)
        history = prefix * "_mission.json"
        write(history, "previous mission history: keep exactly\n")
        history_bytes = read(history)
        manifest_path = prefix * ".manifest.toml"
        for fail_manifest in (false, true)
            old_manifest = Dict("schema_version" => 1, "marker" => "previous run",
                "mission_history" => Dict("schema" => 1))
            open(manifest_path, "w") do io
                TOML.print(io, old_manifest)
            end
            write(prefix * ".feather", "previous result")
            table = DataFrame(time=[0.0, 1.0])
            if fail_manifest
                # An unsupported schema value fails TOML serialization after
                # Feather replacement, exercising a partial bundle write.
                @test_throws ErrorException SM.IOOutputs._write_results_bundle!(table, table.time, config, nothing)
                delete!(old_manifest, "mission_history")
                @test TOML.parsefile(manifest_path) == old_manifest
            else
                SM.IOOutputs._write_results_bundle!(table, table.time, config, 1)
            end
            @test isequal(DataFrame(Arrow.Table(prefix * ".feather")), table)
            @test !haskey(TOML.parsefile(manifest_path), "mission_history")
            @test read(history) == history_bytes
            @test !isfile(prefix * ".csv")
        end
    end
end

@testset "Scene-only simulation replacement drops old mission ownership" begin
    mktempdir() do directory
        output = joinpath(directory, "mission")
        plan = MissionPlan("coast", () -> configuration(output; save_csv=false);
            segments=(MissionSegment(:coast),))
        mission = withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
            run_mission(plan; output_dir=output, save_scene=true, save_payload=false)
        end
        manifest_path = mission.prefix * ".manifest.toml"
        @test haskey(TOML.parsefile(manifest_path), "mission_history")
        history = read(mission.history_path)
        feather = read(mission.prefix * ".feather")
        withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
            memory = run_simulation(configuration(output; results=false, save_csv=false, save_scene=true);
                return_results=true)
            @test isempty(memory.files)
            @test haskey(TOML.parsefile(manifest_path), "mission_history")
            no_outputs = run_simulation(configuration(output; save_csv=false, save_scene=false);
                return_results=true)
            @test isempty(no_outputs.files)
            @test haskey(TOML.parsefile(manifest_path), "mission_history")
            scene_only = run_simulation(configuration(output; save_csv=false, save_scene=true);
                return_results=true)
            @test scene_only.files == [mission.prefix * "_scene.json"]
        end
        @test !haskey(TOML.parsefile(manifest_path), "mission_history")
        @test read(mission.history_path) == history
        @test read(mission.prefix * ".feather") == feather
        @test !isfile(mission.prefix * ".csv")
    end
end

end # module OutputPolicyTests
