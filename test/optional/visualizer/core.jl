using Test, Pkg, SpaceAGORA
include("fixtures.jl")
const VS = SpaceAGORA.VisualizerServices
@testset "Core installation without the renderer" begin
    @test Base.find_package("SpaceAGORAVisualizer") === nothing
    @test all(p.name != "SpaceAGORAVisualizer" for p in values(Pkg.dependencies()))
    @test all(id.name != "SpaceAGORAVisualizer" for id in keys(Base.loaded_modules))
    @test VS.CONTRACT_VERSION == 1
    @test !VS.available() && !visualizer_available()
    @test_throws VisualizerUnavailableError VS.require_available()
    message = try VS.require_available(); "" catch error; sprint(showerror,error) end
    @test occursin("SpaceAGORAVisualizer",message)
    @test occursin("using",message)
    mktempdir() do dir
        f=saved_visualizer_fixture(joinpath(dir,"saved"))
        json=export_viewer_payload(f.prefix; textures=false)
        @test json == f.prefix * "_viewer.json"
        payload=JSON.parsefile(json)
        @test payload["viewer_schema"] == 1
        @test payload["scene"]["schema"] == 1
        @test payload["frames"]["count"] == 3
        @test payload["frames"]["sats"] == 1
        @test !isfile(f.prefix * "_viewer.html")
        for protected in (f.prefix * "_scene.json", f.prefix * ".feather")
            bytes=read(protected)
            @test_throws ArgumentError export_viewer_payload(f.prefix;out=protected,textures=false)
            @test read(protected) == bytes
        end
        @test_throws ArgumentError export_viewer_payload(f.prefix;out=dir,textures=false)
        preserved=joinpath(dir,"existing.html");write(preserved,"keep")
        @test_throws VisualizerUnavailableError export_visualization(f.prefix;out=preserved)
        @test read(preserved,String) == "keep"
        @test_throws VisualizerUnavailableError export_ensemble_visualization(dir;out=preserved)
        @test read(preserved,String) == "keep"
        @test_throws VisualizerUnavailableError run_cli(["visualize","--run=$(f.prefix)","--out=$preserved"])
        @test read(preserved,String) == "keep"

        probe = joinpath(dir, "must-not-run.jl")
        marker = joinpath(dir, "child-started")
        write(probe, "write($(repr(marker)), \"unexpected child\")")
        cli_output = joinpath(dir, "refused-cli")
        withenv("SPACEAGORA_VISUALIZER_ENV" => nothing) do
            @test_throws ArgumentError run_cli(["run", "--example=$probe", "--visualize", "--output-dir=$cli_output"])
            @test_throws ArgumentError run_cli(["run", "--example=$probe", "--visualize", "--print-only"])
        end
        @test !ispath(cli_output)
        @test !ispath(marker)

        destination=joinpath(dir,"refused")
        args=visualizer_args(destination)
        original=deepcopy(args.dynamics_model.spacecraft[1].initial_condition)
        @test_throws VisualizerUnavailableError run_simulation(args;visualization=true)
        @test !ispath(destination)
        @test !args.simulation_settings.save_visualization_scene
        @test isequal(args.dynamics_model.spacecraft[1].initial_condition,original)
        withenv("SPACEAGORA_VISUALIZATION"=>"1") do
            @test_throws VisualizerUnavailableError run_simulation(args)
        end
        @test !ispath(destination)

        campaign=joinpath(dir,"refused-campaign")
        builds=Ref(0);consumed=Ref(0)
        builder=seed -> (builds[]+=1; visualizer_args(destination))
        seeds=((consumed[]+=1;seed) for seed in 1:1)
        @test_throws VisualizerUnavailableError run_monte_carlo_visualization(builder,seeds,campaign)
        @test builds[] == consumed[] == 0
        @test !ispath(campaign)

        # Producing results and scene data remains available without rendering.
        args=with_visualization_scene(visualizer_args(joinpath(dir,"core-run")),true)
        solution=run_simulation(args;visualization=false,return_solution=true)
        @test solution.t[end] ≈ 0.2
        prefix=joinpath(args.simulation_settings.results_directory,"simulation_results")
        @test isfile(prefix * "_scene.json")
        @test isfile(export_viewer_payload(prefix;textures=false))
        @test !isfile(prefix * "_viewer.html")
        result,page=run_monte_carlo_visualization(seed->visualizer_args(destination),[1],
            joinpath(dir,"data-campaign");export_page=false,threads=1)
        @test length(result.successful) == 1
        @test page === nothing
        @test !VS.available()
    end
end
