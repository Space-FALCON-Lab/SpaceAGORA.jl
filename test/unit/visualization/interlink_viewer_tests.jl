module InterlinkViewerTests

using Test
using Arrow
using Base64
using DataFrames
using JSON
using SpaceAGORA
include(joinpath(@__DIR__, "..", "..", "..", "examples", "gve_sma_interlinks.jl"))

const SV = SpaceAGORA.SimulationModel.SceneVisualization

function with_settings(args; overrides...)
    settings = args.simulation_settings
    names = fieldnames(typeof(settings))
    fields = NamedTuple{names}(map(name -> getfield(settings, name), names))
    return SM.SimConfig._with_configuration(args;
        simulation_settings=SimulationSettings(; merge(fields, (; overrides...))...))
end
decode(block) = collect(reinterpret(Float64, base64decode(block)))

# A schedule of two links over three spacecraft, with contact rows.
function hand_schedule()
    spacecraft = [(n_terminal=2,) for _ in 1:3]
    model = InterlinkModel()
    laser = register_link!(model, spacecraft, (1, 1), (2, 1),
        LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=100.0))
    data = register_link!(model, spacecraft, (2, 2), (3, 1),
        DataLink(transmit_power_w=1.0, budget=FixedRateBudget(1e6), range_m=1e6))
    push!(model.history,
        InterLinkHistory(0.0, [laser], [laser, data], [1.0]),
        InterLinkHistory(10.0, [laser, data], [laser, data], [2.0, 5.0]),
        InterLinkHistory(20.0, [data], [data], [6.0]),
        InterLinkHistory(30.0, InterLinkKey[], [data], Float64[]),
        InterLinkHistory(40.0, [laser], InterLinkKey[], [3.0]))
    return link_schedule(model)
end

@testset "Interlinks in the viewer" begin
    @testset "the payload block" begin
        schedule = hand_schedule()
        block = SV.interlink_payload(schedule, 3)
        @test block["schema"] == 1
        @test block["history_start_s"] == 0.0 && block["history_end_s"] == 40.0
        @test [link["link"] for link in block["links"]] == [1, 2]
        @test block["links"][1] == Dict("link" => 1, "kind" => "LaserMomentumLink", "directed" => false,
            "a" => 1, "a_terminal" => 1, "b" => 2, "b_terminal" => 1)
        @test block["links"][2]["kind"] == "DataLink{FixedRateBudget}" && block["links"][2]["directed"]
        @test (block["links"][2]["a"], block["links"][2]["a_terminal"], block["links"][2]["b"]) == (2, 2, 3)
        active = block["active"]
        @test active["count"] == 3
        @test active["index"] == [0, 1, 0]
        @test decode(active["start_s"]) == [0.0, 10.0, 40.0]
        @test decode(active["end_s"]) == [20.0, 30.0, 40.0]
        contact = block["contact"]
        @test contact["count"] == 2 && contact["index"] == [0, 1]
        @test decode(contact["start_s"]) == [0.0, 0.0] && decode(contact["end_s"]) == [20.0, 40.0]
        # Filtering retains recording metadata: no contacts in a recorded span
        # differs from a schedule that never recorded contacts.
        active_only = schedule[schedule.window .== "active", :]
        @test SV.interlink_payload(active_only, 3)["contact"]["count"] == 0
        DataFrames.metadata!(active_only, "record_available", false; style=:note)
        @test SV.interlink_payload(active_only, 3)["contact"] === nothing
        # An empty schedule is an empty block.
        empty_block = SV.interlink_payload(link_schedule(InterlinkModel()), 3)
        @test isempty(empty_block["links"]) && empty_block["active"]["count"] == 0
        recorded_empty = link_schedule(InterlinkModel(record_available=true))
        recorded_block = SV.interlink_payload(recorded_empty, 3)
        @test recorded_block["contact"]["count"] == 0
        @test recorded_block["history_start_s"] === nothing && recorded_block["history_end_s"] === nothing
        # Spacecraft must exist in the scene, and the columns must be there.
        @test_throws ArgumentError SV.interlink_payload(schedule, 2)
        @test_throws ArgumentError SV.interlink_payload(select(schedule, Not(:end_s)), 3)
    end

    @testset "long and streamed viewer input association" begin
        for checkpoints in (false, true)
            mktempdir() do directory
                args = with_settings(gve_sma_configuration(duration=30.0); results=true,
                    results_directory=directory, results_layout=:long, save_csv=false,
                    save_visualization_scene=true, checkpoint_enabled=checkpoints,
                    checkpoint_interval_s=10.0, checkpoint_directory=directory)
                result = run_simulation(args; return_results=true)
                prefix = joinpath(directory, "simulation_results")
                scene = SV.read_visualization_scene(prefix * "_scene.json")
                @test scene.results_layout === :long
                @test scene.link_schedule == "simulation_results_link_schedule.feather"
                @test SV._scene_from_dict(SV.scene_dict(scene)) == scene
                payload = JSON.parsefile(SV.export_viewer_payload(prefix; textures=false))
                @test payload["interlinks"] == JSON.parse(JSON.json(SV.interlink_payload(link_schedule(result), 3)))
                @test isfile(SV.write_viewer_dev_payload(prefix, prefix * "_valid.js"; textures=false))
                trajectory_path = joinpath(directory, scene.results_feather)
                original = read(trajectory_path)
                rows = DataFrame(Arrow.Table(trajectory_path); copycols=true)
                rows.pos_1 .+= 1000.0
                Arrow.write(trajectory_path, rows)
                @test_throws ArgumentError SV.export_viewer_payload(prefix; textures=false, out=prefix * "_stale.json")
                @test_throws ArgumentError SV.write_viewer_dev_payload(prefix, prefix * "_stale.js"; textures=false)
                @test !isfile(prefix * "_stale.json") && !isfile(prefix * "_stale.js")
                @test !haskey(JSON.parsefile(SV.export_viewer_payload(prefix; textures=false, interlinks=false,
                    out=prefix * "_optout.json")), "interlinks")
                @test isfile(SV.write_viewer_dev_payload(prefix, prefix * "_optout.js"; textures=false, interlinks=false))
                write(trajectory_path, original)
                values = ntuple(i -> getfield(scene, i), 10)
                @test SV.VisualizationScene(values...).results_layout === :wide
                @test SV.VisualizationScene(values...).link_schedule === nothing
                @test SV.VisualizationScene(values..., :long).results_layout === :long
                @test SV.VisualizationScene(values..., :long).link_schedule === nothing
                @test SV.VisualizationScene(values..., "schedule.feather").link_schedule == "schedule.feather"
                @test SV.VisualizationScene(values..., "schedule.feather").results_layout === :wide
            end
        end
    end

    @testset "scene, export and checks" begin
        mktempdir() do directory
            args = with_settings(gve_sma_configuration(duration=600.0); results=true, results_directory=directory,
                save_csv=true, save_visualization_scene=true)
            result = run_simulation(args; return_results=true)
            prefix = joinpath(directory, "simulation_results")
            scene_json = JSON.parsefile(prefix * "_scene.json")
            @test scene_json["results"]["link_schedule"] == Dict("file" => "simulation_results_link_schedule.feather")
            scene = SV.read_visualization_scene(prefix * "_scene.json")
            @test scene.link_schedule == "simulation_results_link_schedule.feather"
            @test SV.read_visualization_scene(prefix * "_scene.json") == scene

            path = SV.export_viewer_payload(prefix; textures=false)
            payload = JSON.parsefile(path)
            schedule = link_schedule(result)
            block = payload["interlinks"]
            @test block == JSON.parse(JSON.json(SV.interlink_payload(schedule, 3)))
            @test block["active"]["count"] == count(==("active"), schedule.window)
            @test decode(block["active"]["start_s"]) == schedule.start_s[schedule.window .== "active"]
            @test payload["viewer_schema"] == 1 && payload["scene"]["schema"] == 1

            # The links can be left out; the rest of the payload is the same.
            without = JSON.parsefile(SV.export_viewer_payload(prefix; textures=false, interlinks=false,
                out=prefix * "_without_links.json"))
            @test !haskey(without, "interlinks")
            @test without["frames"] == payload["frames"]

            # The links share the data budget, and the export refuses a budget they exceed.
            @test_throws ArgumentError SV.export_viewer_payload(prefix; textures=false, data_budget_mb=1.0e-6,
                out=prefix * "_tiny.json")
            # The links' bytes come out of the frames' budget: with the links, the frames
            # are those of a payload without them and a budget smaller by their bytes. A
            # long schedule (this one, repeated) makes that difference many frames.
            rows = DataFrame(Arrow.Table(prefix * ".feather"))
            long_schedule = repeat(schedule, 4000)
            link_mb = SV._interlink_payload_bytes(SV.interlink_payload(long_schedule, 3)) / 1.0e6
            budget = link_mb + 0.003
            linked = SV.viewer_payload(scene, rows; link_schedule=long_schedule, data_budget_mb=budget,
                include_textures=false)
            @test linked["frames"] == SV.viewer_payload(scene, rows; data_budget_mb=budget - link_mb,
                include_textures=false)["frames"]
            @test linked["frames"]["count"] < SV.viewer_payload(scene, rows; data_budget_mb=budget,
                include_textures=false)["frames"]["count"]
            @test linked["interlinks"]["active"]["count"] == 4000 * count(==("active"), schedule.window)

            # A schedule that no longer matches its manifest is refused.
            schedule_path = joinpath(directory, "simulation_results_link_schedule.feather")
            original = read(schedule_path)
            write(schedule_path, vcat(original, UInt8[0x00]))
            @test_throws ArgumentError SV.export_viewer_payload(prefix; textures=false, out=prefix * "_stale.json")
            @test haskey(JSON.parsefile(SV.export_viewer_payload(prefix; textures=false, interlinks=false,
                out=prefix * "_stale.json")), "frames")
            write(schedule_path, original)

            # An unchanged schedule cannot be attached to a different valid trajectory.
            trajectory_path = prefix * ".feather"
            trajectory_bytes = read(trajectory_path)
            replacement = DataFrame(Arrow.Table(trajectory_path); copycols=true)
            replacement.sc1_pos_1 .+= 1000.0
            Arrow.write(trajectory_path, replacement)
            @test_throws ArgumentError SV.export_viewer_payload(prefix; textures=false, out=prefix * "_other.json")
            @test_throws ArgumentError SV.write_viewer_dev_payload(prefix, prefix * "_other.js"; textures=false)
            @test !isfile(prefix * "_other.json") && !isfile(prefix * "_other.js")
            @test !haskey(JSON.parsefile(SV.export_viewer_payload(prefix; textures=false, interlinks=false,
                out=prefix * "_without_stale_links.json")), "interlinks")
            @test isfile(SV.write_viewer_dev_payload(prefix, prefix * "_without_stale_links.js"; textures=false, interlinks=false))
            write(trajectory_path, trajectory_bytes)

            # The export never writes over the schedule it reads.
            @test_throws ArgumentError SV.export_viewer_payload(prefix; textures=false, out=schedule_path)
        end
        # A run without interlinks: no schedule in the scene, no block in the payload.
        mktempdir() do directory
            args = with_settings(gve_sma_configuration(duration=60.0, with_interlinks=false); results=true,
                results_directory=directory, save_csv=false, save_visualization_scene=true)
            run_simulation(args)
            prefix = joinpath(directory, "simulation_results")
            @test !haskey(JSON.parsefile(prefix * "_scene.json")["results"], "link_schedule")
            @test SV.read_visualization_scene(prefix * "_scene.json").link_schedule === nothing
            @test !haskey(JSON.parsefile(SV.export_viewer_payload(prefix; textures=false)), "interlinks")
        end
    end
end

end # module InterlinkViewerTests
