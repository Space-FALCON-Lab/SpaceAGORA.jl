using Test, Pkg, TOML
length(ARGS) == 1 || error("Pass core-first or visualizer-first")
if only(ARGS) == "core-first"
    using SpaceAGORA
    @test !SpaceAGORA.VisualizerServices.available()
    using SpaceAGORAVisualizer
elseif only(ARGS) == "visualizer-first"
    using SpaceAGORAVisualizer
    @test all(id.name != "SpaceAGORA" for id in keys(Base.loaded_modules))
    using SpaceAGORA
else
    error("Pass core-first or visualizer-first")
end
include("fixtures.jl")
const VS=SpaceAGORA.VisualizerServices
const V=SpaceAGORAVisualizer
@testset "Explicit renderer integration: $(only(ARGS))" begin
    @test Base.get_extension(SpaceAGORA,:SpaceAGORAVisualizerExt) !== nothing
    @test Base.PkgId(V).uuid == Base.UUID("a4677e92-0b7f-4d10-9ebd-8cde6a49f03a")
    dependency=Pkg.dependencies()[Base.PkgId(V).uuid]
    @test dependency.version == v"0.2.0"
    override=get(ENV,"SPACEAGORA_VISUALIZER_PATH","")
    if isempty(override)
        pin=TOML.parsefile(joinpath(@__DIR__,"..","..","..","packages","VisualizerSource.toml"))
        @test dependency.is_tracking_repo
        @test dependency.git_revision == pin["rev"]
        @test dependency.git_source == pin["url"]
    else
        @test dependency.is_tracking_path
        @test samefile(dependency.source,override)
    end
    @test VS.CONTRACT_VERSION == V.CONTRACT_VERSION == 1
    @test VS.available() && visualizer_available()
    @test VS.require_available() === Base.get_extension(SpaceAGORA,:SpaceAGORAVisualizerExt)
    @test samefile(VS.viewer_directory(),V.viewer_directory())
    mktempdir() do dir
        f=saved_visualizer_fixture(joinpath(dir,"saved"))
        json=export_viewer_payload(f.prefix;textures=false,channels=[(column="metric",)])
        payload=JSON.parsefile(json)
        original=deepcopy(payload)
        @test V.validate_payload(payload) === payload
        standalone=V.visualize(json;out=joinpath(dir,"standalone.html"))
        @test embedded_visualizer_payload(standalone) == payload
        @test payload == original
        page=export_visualization(f.prefix;out=joinpath(dir,"core.html"),textures=false,
            channels=[(column="metric",)],title="<Core> integration")
        rendered=embedded_visualizer_payload(page)
        @test rendered["frames"] == payload["frames"]
        @test rendered["scene"] == payload["scene"]
        @test occursin("&lt;Core&gt; integration",read(page,String))
        @test collect(reinterpret(Float64,base64decode(payload["frames"]["t_s"]))) == [0.0,1.0,2.0]
        @test only(payload["frames"]["channels"])["name"] == "metric"
        @test collect(reinterpret(Float32,base64decode(only(payload["frames"]["channels"])["data"]))) == Float32[3,2,1]
        cli=joinpath(dir,"cli.html")
        @test run_cli(["visualize","--run=$(f.prefix)","--no-textures","--out=$cli"];io=IOBuffer()) == 0
        @test embedded_visualizer_payload(cli)["scene"] == payload["scene"]

        # The actual CLI child must retain the explicit renderer even when a
        # normal example continues in the examples environment. Julia -L must
        # preserve its entrypoint.
        project = dirname(Base.active_project())
        CLI = SpaceAGORA.SpaceAGORACLI
        probe_dir = joinpath(dir, "child path with spaces")
        mkpath(probe_dir)
        probe = joinpath(probe_dir, "viewer probe.jl")
        write(probe, """
            using SpaceAGORA
            @assert visualizer_available()
            @assert abspath(PROGRAM_FILE) == @__FILE__
            @assert isempty(ARGS)
            include($(repr(joinpath(CLI.REPO_ROOT, "examples", "common.jl"))))
            @assert Base.active_project() == $(repr(joinpath(CLI.REPO_ROOT, "examples", "Project.toml")))
            @assert visualizer_available()
            module ChildFixture
                include($(repr(joinpath(@__DIR__, "fixtures.jl"))))
            end
            run_simulation(ChildFixture.visualizer_args(ENV["SPACEAGORA_CLI_OUTPUT_DIR"]))
            println("VISUALIZER_CHILD_OK")
            """)
        child_output = joinpath(dir, "cli-child-run")
        arguments = ["run", "--example=$probe", "--output-dir=$child_output", "--visualize"]
        withenv("SPACEAGORA_VISUALIZER_ENV" => project) do
            @test CLI._visualizer_project() == project
            printed = sprint(io -> @test run_cli(vcat(arguments, ["--print-only"]); io, errio=io) == 0)
            @test occursin("project=$project\n", printed)
            @test occursin("visualizer_preload.jl", printed)
            @test !ispath(child_output)
            log = joinpath(dir, "cli-child.log")
            status = open(log, "w") do io
                run_cli(arguments; io, errio=io)
            end
            @test status == 0
            @test occursin("VISUALIZER_CHILD_OK", read(log, String))
            child_page = joinpath(child_output, "simulation_results_viewer.html")
            @test embedded_visualizer_payload(child_page)["frames"]["count"] > 0
        end
        withenv("SPACEAGORA_VISUALIZER_ENV" => nothing) do
            @test CLI._visualizer_project() == project
        end
        invalid_env = joinpath(dir, "invalid-environment")
        mkpath(invalid_env)
        write(joinpath(invalid_env, "Project.toml"), "[deps]\nSpaceAGORA = \"afbfb69f-5c0b-4832-b760-43725dff8540\"\nSpaceAGORAVisualizer = \"00000000-0000-0000-0000-000000000000\"\n")
        withenv("SPACEAGORA_VISUALIZER_ENV" => invalid_env) do
            @test_throws ArgumentError CLI._visualizer_project()
        end
        child_argv = collect(CLI._child_command(probe, ["first argument", "--second=two words"];
            project, visualize=true).exec)
        @test child_argv[end-2:end] == [probe, "first argument", "--second=two words"]
        @test count(==("--project=$project"), child_argv) == 1

        # The direct script entrypoint explicitly loads the prepared companion.
        entry_page = joinpath(dir, "entrypoint.html")
        entry = joinpath(CLI.REPO_ROOT, "src", "interfaces", "cli", "main.jl")
        command = `$(Base.julia_cmd()) --startup-file=no --project=$project $entry visualize --run=$(f.prefix) --no-textures --out=$entry_page`
        entry_log = joinpath(dir, "entrypoint.log")
        entry_status = withenv("SPACEAGORA_VISUALIZER_ENV" => nothing) do
            open(entry_log, "w") do io
                run(pipeline(ignorestatus(command); stdout=io, stderr=io)).exitcode
            end
        end
        @test entry_status == 0
        @test embedded_visualizer_payload(entry_page)["scene"] == payload["scene"]
        invalid=deepcopy(payload);invalid["viewer_schema"]=2
        @test_throws V.CompatibilityError V.render_viewer_html(invalid)
        @test payload == original
        args=visualizer_args(joinpath(dir,"run"))
        run_simulation(args;visualization=true)
        @test !args.simulation_settings.save_visualization_scene
        @test isfile(joinpath(args.simulation_settings.results_directory,"simulation_results_viewer.html"))
        result,ensemble=run_monte_carlo_visualization(seed->visualizer_args(joinpath(dir,"unused")),[1],
            joinpath(dir,"campaign");threads=1,textures=false)
        @test length(result.successful) == 1
        @test embedded_visualizer_payload(ensemble)["ensemble"]["count"] == 1

        builder_calls = Ref(0)
        mission = MissionPlan("Saved viewer mission", () -> begin
            builder_calls[] += 1
            visualizer_args(joinpath(dir, "unused-mission-builder"))
        end)
        mission_run = run_mission(mission; output_dir=joinpath(dir, "mission"))
        @test builder_calls[] == 1
        mission_payload = JSON.parsefile(mission_run.payload_path)
        mission_page = V.visualize(mission_run)
        @test embedded_visualizer_payload(mission_page) == mission_payload
        @test builder_calls[] == 1
        moved = joinpath(dir, "moved-mission")
        mv(mission_run.directory, moved)
        @test !ispath(mission_run.directory)
        moved_page = V.visualize(moved)
        @test embedded_visualizer_payload(moved_page) == mission_payload
        @test builder_calls[] == 1
    end
end
