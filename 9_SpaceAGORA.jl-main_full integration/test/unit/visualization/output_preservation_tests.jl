module ViewerOutputPreservationTests

using Test
using TOML
using OrdinaryDiffEq: DiscreteCallback
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "helpers", "load_visualizer.jl"))
include(joinpath(@__DIR__, "..", "..", "optional_visualizer", "fixtures.jl"))

# Most cases use tiny saved tables; the overwrite regression also runs the real
# mission and simulation entry points against the same output directory.
function bind_fixture_history(prefix, history)
    write(prefix * "_mission.json", JSON.json(history))
    manifest = Dict("schema_version" => 1, "files" => Dict("feather" => Dict(
        "sha256" => SM.IOSerialization._sha256_hex(prefix * ".feather"))))
    open(prefix * ".manifest.toml", "w") do io
        TOML.print(io, manifest)
    end
    SV._bind_mission_history!(prefix)
end

function development_payload(prefix, output; kwargs...)
    path = write_viewer_dev_payload(prefix, output; textures=false, kwargs...)
    return JSON.parse(chop(replace(read(path, String), "window.SPACEAGORA_VIEWER = " => ""); tail=2))
end
function aliases(path, directory, label)
    symbolic = joinpath(directory, label * "-symbolic")
    linked = joinpath(directory, label * "-hard")
    symlink(path, symbolic)
    hardlink(path, linked)
    return (path, symbolic, linked)
end

function preserves_inputs(writer, target, originals)
    @test_throws ArgumentError writer(target)
    for (path, bytes) in originals
        @test isfile(path)
        @test read(path) == bytes
    end
end

@testset "Single-run exporters preserve saved inputs and directories" begin
    mktempdir() do directory
        saved = saved_visualizer_fixture(joinpath(directory, "saved"))
        history = Dict("schema" => 1, "name" => "coast", "segments" => [Dict("name" => "coast", "start_s" => 0.0, "end_s" => 2.0)], "transitions" => [])
        bind_fixture_history(saved.prefix, history)
        payload_path = export_viewer_payload(saved.prefix; textures=false)
        @test JSON.parsefile(payload_path)["mission"] == history
        inputs = (saved.prefix * "_scene.json", saved.prefix * ".feather", saved.prefix * "_mission.json", saved.prefix * ".manifest.toml")
        originals = Dict(path => read(path) for path in (inputs..., payload_path))
        sentinel = joinpath(directory, "keep.txt")
        write(sentinel, "keep exactly")
        directory_alias = joinpath(directory, "saved-directory-alias")
        symlink(dirname(saved.prefix), directory_alias)
        writers = (
            payload = out -> export_viewer_payload(saved.prefix; out, textures=false),
            html = out -> export_visualization(saved.prefix; out, textures=false),
            development = out -> write_viewer_dev_payload(saved.prefix, out; textures=false),
        )
        targets = [aliases(path, directory, "input-$i") for (i,path) in enumerate(inputs)]
        for writer in values(writers)
            for target in (directory, directory_alias)
                preserves_inputs(writer, target, originals)
                @test isdir(target)
                @test read(sentinel, String) == "keep exactly"
            end
            for group in targets, target in group
                preserves_inputs(writer, target, originals)
            end
        end
        for target in aliases(payload_path, directory, "portable")
            preserves_inputs(writers.html, target, originals)
            preserves_inputs(writers.development, target, originals)
        end
        @test islink(directory_alias)
        @test all(islink(group[2]) && samefile(group[1], group[3]) for group in targets)

        # Explicit replacement of separate output files remains useful.
        json = joinpath(directory, "separate.json")
        html = joinpath(directory, "separate.html")
        javascript = joinpath(directory, "dev_data.js")
        for target in (json, html, javascript)
            write(target, "replace this output")
        end
        @test writers.payload(json) == json
        @test JSON.parsefile(json) == JSON.parsefile(payload_path)
        @test writers.html(html) == html
        @test embedded_visualizer_payload(html) == JSON.parsefile(payload_path)
        @test writers.development(javascript) == javascript
        @test startswith(read(javascript, String), "window.SPACEAGORA_VIEWER = ")
        @test read(sentinel, String) == "keep exactly"
    end
end

@testset "Ensemble HTML protects every sample and its manifest" begin
    mktempdir() do directory
        campaign = joinpath(directory, "campaign")
        fixtures = [saved_visualizer_fixture(joinpath(campaign, "sample_$(lpad(i,4,'0'))")) for i in 1:2]
        samples = [EnsembleSample(i, string(i), true, Float64(i), "sample $i", dirname(f.prefix)) for (i,f) in enumerate(fixtures)]
        manifest = write_ensemble_manifest(campaign, samples; scalar_name="synthetic fixture")
        inputs = [manifest]
        for fixture in fixtures
            history = Dict("schema" => 1, "name" => "coast", "segments" => [Dict("name" => "coast", "start_s" => 0.0, "end_s" => 2.0)], "transitions" => [])
            bind_fixture_history(fixture.prefix, history)
            push!(inputs, fixture.prefix * "_scene.json", fixture.prefix * ".feather", fixture.prefix * "_mission.json", fixture.prefix * ".manifest.toml",
                export_viewer_payload(fixture.prefix; textures=false))
        end
        originals = Dict(path => read(path) for path in inputs)
        sentinel = joinpath(campaign, "keep.txt"); write(sentinel, "keep exactly")
        directory_alias = joinpath(directory, "campaign-alias"); symlink(campaign, directory_alias)
        writer = out -> export_ensemble_visualization(campaign; out, textures=false, max_frames=3)
        for target in (campaign, directory_alias)
            preserves_inputs(writer, target, originals)
            @test isdir(target)
            @test read(sentinel, String) == "keep exactly"
        end
        for (i,path) in enumerate(inputs), target in aliases(path, directory, "ensemble-$i")
            preserves_inputs(writer, target, originals)
        end
        @test islink(directory_alias)
        output = joinpath(directory, "ensemble.html"); write(output, "replace this output")
        @test writer(output) == output
        @test embedded_visualizer_payload(output)["ensemble"]["count"] == 2
        @test read(sentinel, String) == "keep exactly"
    end
end


@testset "Mission history export and boundary sampling" begin
    mktempdir() do directory
        saved = saved_visualizer_fixture(joinpath(directory, "saved"))
        table = DataFrame(time=collect(0.0:100.0), sc1_pos_1=collect(0.0:100.0) .* 1000,
            sc1_pos_2=zeros(101), sc1_pos_3=zeros(101))
        mission = Dict{String, Any}("schema" => 1, "name" => "Coast and maneuver",
            "segments" => [Dict("name" => "coast", "start_s" => 0.0, "end_s" => 37.0),
                Dict("name" => "maneuver", "start_s" => 37.0, "end_s" => 100.0)],
            "transitions" => [Dict("time_s" => 37.0, "from" => "coast", "to" => "maneuver", "reason" => "deadline")])
        before = deepcopy(mission)
        payload = SV.viewer_payload(saved.scene, table; mission, include_textures=false, max_frames=8)
        times = collect(reinterpret(Float64, base64decode(payload["frames"]["t_s"])))
        positions = collect(reinterpret(Float64, base64decode(payload["frames"]["pos_km"])))
        @test payload["mission"] == mission == before
        @test payload["mission"] !== mission
        @test payload["frames"]["count"] <= 8
        @test all(t -> t in times, [0.0, 37.0, 100.0])
        @test positions[1:3:end] == times
        @test !haskey(SV.viewer_payload(saved.scene, table; include_textures=false), "mission")
        @test_throws ArgumentError SV.viewer_payload(saved.scene, table; mission, include_textures=false, max_frames=2)
        @test_throws ArgumentError SV.viewer_payload(saved.scene, table; mission, include_textures=false, data_budget_mb=0.000064)
        three = SV.viewer_payload(saved.scene, table; mission, include_textures=false, max_frames=3)
        @test collect(reinterpret(Float64, base64decode(three["frames"]["t_s"]))) == [0.0, 37.0, 100.0]
        four = SV.viewer_payload(saved.scene, table; mission, include_textures=false, max_frames=4)
        @test four["frames"]["count"] == 4
        absent = deepcopy(mission)
        absent["segments"][1]["end_s"] = absent["segments"][2]["start_s"] = absent["transitions"][1]["time_s"] = 37.5
        without_row = SV.viewer_payload(saved.scene, table; mission=absent, include_textures=false, max_frames=8)
        @test without_row["mission"] == absent
        @test !(37.5 in reinterpret(Float64, base64decode(without_row["frames"]["t_s"])))
        duplicate = deepcopy(mission); duplicate["segments"][2]["name"] = "coast"
        @test_throws ArgumentError SV.viewer_payload(saved.scene, table; mission=duplicate, include_textures=false)
        bad_transition = deepcopy(mission); bad_transition["transitions"][1]["time_s"] = 38.0
        @test_throws ArgumentError SV.viewer_payload(saved.scene, table; mission=bad_transition, include_textures=false)
        for change in (
            m -> m["schema"] = true,
            m -> m["schema"] = 2,
            m -> m["name"] = " ",
            m -> m["segments"] = "coast",
            m -> m["segments"] = [nothing],
            m -> m["segments"][1]["start_s"] = -1,
            m -> m["segments"][1]["end_s"] = Inf,
            m -> m["segments"][1]["end_s"] = true,
            m -> m["segments"][2]["start_s"] = 38,
            m -> m["transitions"] = "deadline",
            m -> m["transitions"] = [nothing],
            m -> empty!(m["transitions"]),
            m -> m["transitions"][1]["from"] = "unknown",
            m -> m["transitions"][1]["to"] = nothing,
            m -> m["transitions"][1]["reason"] = " ",
        )
            invalid = deepcopy(mission); change(invalid)
            @test_throws ArgumentError SV.viewer_payload(saved.scene, table; mission=invalid, include_textures=false)
        end
        @test_throws ArgumentError SV.viewer_payload(saved.scene, table; mission=[], include_textures=false)
        empty_history = deepcopy(mission); empty_history["segments"] = []; empty_history["transitions"] = []
        @test SV.viewer_payload(saved.scene, table; mission=empty_history, include_textures=false)["mission"] == empty_history
        @test SV.viewer_payload(saved.scene, table[1:0, :]; mission, include_textures=false)["frames"]["count"] == 0
        zero_history = deepcopy(mission)
        zero_history["segments"][1]["end_s"] = zero_history["segments"][2]["start_s"] = zero_history["transitions"][1]["time_s"] = 0.0
        @test SV.viewer_payload(saved.scene, table; mission=zero_history, include_textures=false)["mission"] == zero_history
        final_history = deepcopy(mission)
        push!(final_history["transitions"], Dict("time_s" => 100.0, "from" => "maneuver", "to" => nothing, "reason" => "horizon"))
        @test SV.viewer_payload(saved.scene, table; mission=final_history, include_textures=false)["mission"] == final_history
        # The exporter discovers the saved record after the directory is moved.
        Arrow.write(saved.prefix * ".feather", table)
        bind_fixture_history(saved.prefix, mission)
        moved = joinpath(directory, "moved"); mv(dirname(saved.prefix), moved)
        prefix = joinpath(moved, "simulation_results")
        json = export_viewer_payload(prefix; textures=false, max_frames=8)
        @test JSON.parsefile(json)["mission"] == mission
        html = export_visualization(prefix; textures=false, max_frames=8)
        @test embedded_visualizer_payload(html)["mission"] == mission
        js = write_viewer_dev_payload(prefix, joinpath(moved, "dev.js"); textures=false, max_frames=8)
        source = read(js, String)
        @test JSON.parse(chop(replace(source, "window.SPACEAGORA_VIEWER = " => ""); tail=2))["mission"] == mission
        write(prefix * "_mission.json", JSON.json(bad_transition))
        original = read(json)
        @test_throws ArgumentError export_viewer_payload(prefix; textures=false)
        @test read(json) == original
    end
end

@testset "Saved history belongs to one result write" begin
    mktempdir() do directory
        for bundle in (true, false), scene in (true, false)
            output = joinpath(directory, "bundle-$bundle-scene-$scene")
            plan = MissionPlan("recorded coast", () -> visualizer_args(output);
                segments=(MissionSegment(:coast),))
            mission = withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
                run_mission(plan; output_dir=output)
            end
            history_bytes = read(mission.history_path)
            history = JSON.parse(String(copy(history_bytes)))
            @test JSON.parsefile(mission.payload_path)["mission"] == history
            @test development_payload(mission.prefix, joinpath(output, "dev.js"))["mission"] == history
            manifest_path = mission.prefix * ".manifest.toml"
            manifest_bytes = read(manifest_path)
            settings = mission.results.configuration.simulation_settings
            names = fieldnames(typeof(settings))
            original = NamedTuple{names}(map(n -> getfield(settings, n), names))
            function replay(; results=true)
                configured = SM.SimulationSettings(; merge(original,
                    (; results, save_visualization_scene=scene))...)
                args = SM.SimConfig._with_configuration(mission.results.configuration;
                    simulation_settings=configured)
                return withenv("SPACEAGORA_SAVE_BUNDLE" => bundle ? "1" : "0") do
                    run_simulation(args; return_results=true)
                end
            end
            # Returning an in-memory table does not replace a saved run.
            replay(; results=false)
            @test read(manifest_path) == manifest_bytes
            @test read(mission.history_path) == history_bytes
            @test JSON.parsefile(export_viewer_payload(mission.prefix; textures=false))["mission"] == history
            # A normal saved simulation has no operational history, even if it
            # reproduces the same coast and leaves an old scene or Arrow file.
            replay()
            @test !haskey(TOML.parsefile(manifest_path), "mission_history")
            @test read(mission.history_path) == history_bytes
            @test !haskey(JSON.parsefile(export_viewer_payload(mission.prefix; textures=false)), "mission")
            @test !haskey(development_payload(mission.prefix, joinpath(output, "dev.js")), "mission")
            @test !haskey(embedded_visualizer_payload(export_visualization(mission.prefix; textures=false)), "mission")
        end
        output = joinpath(directory, "portable")
        plan = MissionPlan("portable coast", () -> visualizer_args(output);
            segments=(MissionSegment(:coast),))
        mission = withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
            run_mission(plan; output_dir=output)
        end
        history = JSON.parsefile(mission.history_path)
        moved = joinpath(directory, "relocated"); mv(output, moved)
        prefix = joinpath(moved, basename(mission.prefix))
        @test JSON.parsefile(export_viewer_payload(prefix; textures=false))["mission"] == history
        @test development_payload(prefix, joinpath(moved, "dev.js"))["mission"] == history
        @test embedded_visualizer_payload(export_visualization(prefix; textures=false))["mission"] == history
    end
end

@testset "Failed and scene-disabled runs do not inherit mission ownership" begin
    mktempdir() do directory
        output = joinpath(directory, "no-scene")
        plan = MissionPlan("coast without a scene", () -> visualizer_args(output);
            segments=(MissionSegment(:coast),))
        mission = withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
            run_mission(plan; output_dir=output, save_scene=false, save_payload=false)
        end
        @test isfile(mission.history_path)
        @test JSON.parsefile(mission.history_path)["name"] == plan.name
        @test !haskey(TOML.parsefile(mission.prefix * ".manifest.toml"), "mission_history")
        @test !isfile(mission.prefix * "_scene.json")
        for bundle in (true, false)
            output = joinpath(directory, "failure-$bundle")
            plan = MissionPlan("coast before failure", () -> visualizer_args(output);
                segments=(MissionSegment(:coast),))
            mission = withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
                run_mission(plan; output_dir=output)
            end
            history_bytes = read(mission.history_path)
            @test haskey(TOML.parsefile(mission.prefix * ".manifest.toml"), "mission_history")
            failed_after_step = Ref(false)
            callback = _ -> DiscreteCallback((u, t, integrator) -> t > 0.0,
                integrator -> begin
                    failed_after_step[] = true
                    error("injected failure after an accepted step")
                end; save_positions=(false, false))
            args = SM.SimConfig._with_configuration(mission.results.configuration;
                solver_config=SolverConfig(solver_mode=:tsit5, parallel=false))
            @test_throws ErrorException withenv("SPACEAGORA_SAVE_BUNDLE" => bundle ? "1" : "0") do
                run_simulation(args; extra_callback_factories=(callback,))
            end
            @test failed_after_step[]
            @test !haskey(TOML.parsefile(mission.prefix * ".manifest.toml"), "mission_history")
            @test read(mission.history_path) == history_bytes
            @test !haskey(JSON.parsefile(export_viewer_payload(mission.prefix; textures=false)), "mission")
            @test !haskey(development_payload(mission.prefix, joinpath(output, "dev.js")), "mission")
        end
    end
end

@testset "Only verified history is discovered; caller input stays explicit" begin
    mktempdir() do directory
        saved = saved_visualizer_fixture(directory)
        history = Dict("schema" => 1, "name" => "coast", "segments" => [Dict("name" => "coast", "start_s" => 0.0, "end_s" => 2.0)], "transitions" => [])
        history_path = saved.prefix * "_mission.json"
        manifest_path = saved.prefix * ".manifest.toml"
        write(history_path, JSON.json(history))
        history_bytes = read(history_path)
        json = export_viewer_payload(saved.prefix; textures=false)
        @test !haskey(JSON.parsefile(json), "mission")
        @test !haskey(development_payload(saved.prefix, joinpath(directory, "dev.js")), "mission")
        @test read(history_path) == history_bytes
        # A malformed unrelated manifest does not turn a standalone scene into
        # an owned mission run or stop ordinary legacy scene export.
        write(manifest_path, "not = valid = toml")
        @test !haskey(JSON.parsefile(export_viewer_payload(saved.prefix; textures=false)), "mission")
        SM.IOOutputs._write_results_csv!(saved.table, saved.args)
        @test read(manifest_path, String) == "not = valid = toml"
        @test read(history_path) == history_bytes
        bind_fixture_history(saved.prefix, history)
        json = export_viewer_payload(saved.prefix; textures=false)
        original_json = read(json)
        for input in (history_path, saved.prefix * "_scene.json", saved.prefix * ".feather")
            original = read(input)
            if endswith(input, ".feather")
                modified = copy(saved.table); modified.sc1_pos_1 .+= 1
                Arrow.write(input, modified)
            else
                write(input, vcat(original, UInt8[0x20]))
            end
            @test_throws ArgumentError export_viewer_payload(saved.prefix; textures=false)
            @test_throws ArgumentError write_viewer_dev_payload(saved.prefix, joinpath(directory, "dev.js"); textures=false)
            @test read(json) == original_json
            # Explicit history is a deliberate caller choice, independent of
            # automatic saved-run association, in both export entry points.
            explicit = export_viewer_payload(saved.prefix; out=joinpath(directory, "explicit.json"), textures=false, mission=history)
            @test JSON.parsefile(explicit)["mission"] == history
            @test development_payload(saved.prefix, joinpath(directory, "explicit.js"); mission=history)["mission"] == history
            write(input, original)
        end
        # Rewriting exactly the same CSV still invalidates ownership. The
        # durable history and all unrelated manifest fields remain intact.
        manifest = TOML.parsefile(manifest_path)
        expected = deepcopy(manifest); delete!(expected, "mission_history")
        arrow_bytes = read(saved.prefix * ".feather")
        SM.IOOutputs._write_results_csv!(saved.table, saved.args)
        @test TOML.parsefile(manifest_path) == expected
        @test read(history_path) == history_bytes
        @test read(saved.prefix * ".feather") == arrow_bytes
        @test !haskey(JSON.parsefile(export_viewer_payload(saved.prefix; textures=false)), "mission")
    end
end

end # module ViewerOutputPreservationTests
